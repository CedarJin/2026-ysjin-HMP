#include <algorithm>
#include <atomic>
#include <cmath>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <vector>

#ifdef _OPENMP
#include <omp.h>
#endif

struct FastqRecord {
    std::string header;
    std::string sequence;
    std::string plus;
    std::string quality;
};

struct PairInfo {
    uint64_t pair_index;
    uint64_t read1_offset;
    uint64_t read2_offset;
    uint32_t seed;
    uint16_t read1_length;
    uint16_t read2_length;
};

struct Options {
    std::vector<std::string> r1;
    std::vector<std::string> r2;
    std::string out_r1;
    std::string out_r2;
    std::string singleton;
    std::string report;
    int threads = 1;
    int seed_length = 5;
    double max_diff_rate = 0.03;
    int min_mean_quality = 20;
    int trim_quality = 3;
    int min_length = 60;
};

static std::string require_value(int &i, int argc, char **argv) {
    if (++i >= argc) throw std::runtime_error("missing value for command-line option");
    return argv[i];
}

static Options parse_args(int argc, char **argv) {
    Options opt;
    for (int i = 1; i < argc; ++i) {
        const std::string arg = argv[i];
        if (arg == "--r1") opt.r1.push_back(require_value(i, argc, argv));
        else if (arg == "--r2") opt.r2.push_back(require_value(i, argc, argv));
        else if (arg == "--out-r1") opt.out_r1 = require_value(i, argc, argv);
        else if (arg == "--out-r2") opt.out_r2 = require_value(i, argc, argv);
        else if (arg == "--singleton") opt.singleton = require_value(i, argc, argv);
        else if (arg == "--report") opt.report = require_value(i, argc, argv);
        else if (arg == "--threads") opt.threads = std::stoi(require_value(i, argc, argv));
        else if (arg == "--min-identical-bases") opt.seed_length = std::stoi(require_value(i, argc, argv));
        else if (arg == "--max-diff-rate") opt.max_diff_rate = std::stod(require_value(i, argc, argv));
        else if (arg == "--min-mean-quality") opt.min_mean_quality = std::stoi(require_value(i, argc, argv));
        else if (arg == "--trim-quality") opt.trim_quality = std::stoi(require_value(i, argc, argv));
        else if (arg == "--min-length") opt.min_length = std::stoi(require_value(i, argc, argv));
        else throw std::runtime_error("unknown option: " + arg);
    }
    if (opt.r1.empty() || opt.r1.size() != opt.r2.size() || opt.out_r1.empty() || opt.out_r2.empty() ||
        opt.singleton.empty() || opt.report.empty()) {
        throw std::runtime_error("matching --r1/--r2 inputs and all output options are required");
    }
    if (opt.seed_length < 1 || opt.seed_length > 8) {
        throw std::runtime_error("--min-identical-bases must be between 1 and 8");
    }
    return opt;
}

static bool read_fastq(std::ifstream &in, FastqRecord &record) {
    if (!std::getline(in, record.header)) return false;
    if (!std::getline(in, record.sequence) || !std::getline(in, record.plus) ||
        !std::getline(in, record.quality)) {
        throw std::runtime_error("truncated FASTQ record after " + record.header);
    }
    auto strip_cr = [](std::string &value) {
        if (!value.empty() && value.back() == '\r') value.pop_back();
    };
    strip_cr(record.header);
    strip_cr(record.sequence);
    strip_cr(record.plus);
    strip_cr(record.quality);
    if (record.header.empty() || record.header[0] != '@' || record.plus.empty() ||
        record.plus[0] != '+' || record.sequence.size() != record.quality.size()) {
        throw std::runtime_error("malformed FASTQ record " + record.header);
    }
    return true;
}

static std::string read_id(const std::string &header) {
    size_t end = header.find_first_of(" \t");
    std::string id = header.substr(1, end == std::string::npos ? end : end - 1);
    if (id.size() > 2 && id[id.size() - 2] == '/' &&
        (id.back() == '1' || id.back() == '2')) {
        id.resize(id.size() - 2);
    }
    return id;
}

static int base_code(char base) {
    switch (base) {
        case 'A': case 'a': return 0;
        case 'C': case 'c': return 1;
        case 'G': case 'g': return 2;
        case 'T': case 't': return 3;
        default: return -1;
    }
}

static bool quality_eligible(const FastqRecord &record, int seed_length, int min_mean_quality) {
    if (record.sequence.size() < static_cast<size_t>(seed_length)) return false;
    for (int i = 0; i < seed_length; ++i) {
        if (base_code(record.sequence[i]) < 0) return false;
    }
    int64_t quality_sum = 0;
    for (char value : record.quality) quality_sum += static_cast<unsigned char>(value) - 33;
    return quality_sum / static_cast<int64_t>(record.quality.size()) >= min_mean_quality;
}

static uint32_t pair_seed(const FastqRecord &r1, const FastqRecord &r2, int seed_length) {
    uint32_t seed = 0;
    for (int i = 0; i < seed_length; ++i) seed = (seed << 2) | base_code(r1.sequence[i]);
    for (int i = 0; i < seed_length; ++i) seed = (seed << 2) | base_code(r2.sequence[i]);
    return seed;
}

static bool matches(const PairInfo &lhs, const PairInfo &rhs,
                    const std::vector<uint8_t> &bases, int seed_length, double max_diff_rate) {
    const size_t len1 = std::min(lhs.read1_length, rhs.read1_length);
    const size_t len2 = std::min(lhs.read2_length, rhs.read2_length);
    const int max_errors = static_cast<int>(std::floor((len1 + len2) * max_diff_rate));
    int errors = 0;
    for (size_t i = seed_length; i < len1; ++i) {
        if (bases[lhs.read1_offset + i] != bases[rhs.read1_offset + i] && ++errors > max_errors) return false;
    }
    for (size_t i = seed_length; i < len2; ++i) {
        if (bases[lhs.read2_offset + i] != bases[rhs.read2_offset + i] && ++errors > max_errors) return false;
    }
    return true;
}

static size_t bwa_trim_length(const std::string &quality, int threshold) {
    int64_t pos = static_cast<int64_t>(quality.size());
    int64_t max_pos = pos;
    int area = 0;
    int max_area = 0;
    while (pos > 0 && area >= 0) {
        area += threshold - (static_cast<unsigned char>(quality[pos - 1]) - 33);
        if (area > max_area) {
            max_area = area;
            max_pos = pos;
        }
        --pos;
    }
    return static_cast<size_t>(pos == 0 ? 1 : max_pos);
}

static void write_fastq(std::ofstream &out, const FastqRecord &record, size_t length) {
    out << record.header << '\n'
        << record.sequence.substr(0, length) << '\n'
        << record.plus << '\n'
        << record.quality.substr(0, length) << '\n';
}

int main(int argc, char **argv) {
    try {
        const Options opt = parse_args(argc, argv);
        std::vector<PairInfo> pairs;
        std::vector<uint8_t> bases;
        uintmax_t input_bytes = 0;
        for (size_t i = 0; i < opt.r1.size(); ++i) {
            input_bytes += std::filesystem::file_size(opt.r1[i]);
            input_bytes += std::filesystem::file_size(opt.r2[i]);
        }
        bases.reserve(static_cast<size_t>(input_bytes / 3));
        pairs.reserve(static_cast<size_t>(input_bytes / 550));

        FastqRecord rec1, rec2;
        uint64_t total_pairs = 0;
        for (size_t file_index = 0; file_index < opt.r1.size(); ++file_index) {
            std::ifstream r1(opt.r1[file_index]), r2(opt.r2[file_index]);
            if (!r1 || !r2) throw std::runtime_error("cannot open input FASTQ files");
            while (true) {
                const bool have1 = read_fastq(r1, rec1);
                const bool have2 = read_fastq(r2, rec2);
                if (!have1 && !have2) break;
                if (have1 != have2) throw std::runtime_error("mate FASTQ files have different record counts");
                if (read_id(rec1.header) != read_id(rec2.header)) {
                    throw std::runtime_error("mate ID mismatch: " + rec1.header + " vs " + rec2.header);
                }
                if (rec1.sequence.size() > UINT16_MAX || rec2.sequence.size() > UINT16_MAX) {
                    throw std::runtime_error("read length exceeds supported FASTQ limit");
                }
                if (quality_eligible(rec1, opt.seed_length, opt.min_mean_quality) &&
                    quality_eligible(rec2, opt.seed_length, opt.min_mean_quality)) {
                    const uint64_t off1 = bases.size();
                    bases.insert(bases.end(), rec1.sequence.begin(), rec1.sequence.end());
                    const uint64_t off2 = bases.size();
                    bases.insert(bases.end(), rec2.sequence.begin(), rec2.sequence.end());
                    pairs.push_back({total_pairs, off1, off2, pair_seed(rec1, rec2, opt.seed_length),
                                     static_cast<uint16_t>(rec1.sequence.size()),
                                     static_cast<uint16_t>(rec2.sequence.size())});
                }
                ++total_pairs;
                if (total_pairs % 10000000 == 0) std::cerr << "Loaded " << total_pairs << " read pairs\n";
            }
        }

        std::cerr << "Sorting " << pairs.size() << " quality-eligible pairs by "
                  << 2 * opt.seed_length << " seed bases\n";
        const uint64_t quality_eligible_pairs = pairs.size();
        std::sort(pairs.begin(), pairs.end(), [](const PairInfo &a, const PairInfo &b) {
            if (a.seed != b.seed) return a.seed < b.seed;
            return a.pair_index < b.pair_index;
        });

        std::vector<std::pair<size_t, size_t>> groups;
        for (size_t begin = 0; begin < pairs.size();) {
            size_t end = begin + 1;
            while (end < pairs.size() && pairs[end].seed == pairs[begin].seed) ++end;
            groups.emplace_back(begin, end);
            begin = end;
        }

        std::vector<uint8_t> duplicate(total_pairs, 0);
        std::atomic<uint64_t> duplicate_pairs{0};
#ifdef _OPENMP
        omp_set_num_threads(opt.threads);
#pragma omp parallel for schedule(dynamic, 1)
#endif
        for (int64_t group_index = 0; group_index < static_cast<int64_t>(groups.size()); ++group_index) {
            const auto [begin, end] = groups[group_index];
            for (size_t i = begin; i < end; ++i) {
                if (duplicate[pairs[i].pair_index]) continue;
                for (size_t j = i + 1; j < end; ++j) {
                    if (duplicate[pairs[j].pair_index]) continue;
                    if (matches(pairs[i], pairs[j], bases, opt.seed_length, opt.max_diff_rate)) {
                        duplicate[pairs[j].pair_index] = 1;
                        ++duplicate_pairs;
                    }
                }
            }
        }

        pairs.clear();
        pairs.shrink_to_fit();
        bases.clear();
        bases.shrink_to_fit();

        std::ofstream out1(opt.out_r1), out2(opt.out_r2), singleton(opt.singleton);
        if (!out1 || !out2 || !singleton) throw std::runtime_error("cannot open output FASTQ files");

        uint64_t pair_index = 0;
        uint64_t kept_pairs = 0;
        uint64_t singleton_reads = 0;
        uint64_t short_pairs = 0;
        uint64_t quality_bases = 0;
        for (size_t file_index = 0; file_index < opt.r1.size(); ++file_index) {
            std::ifstream r1(opt.r1[file_index]), r2(opt.r2[file_index]);
            if (!r1 || !r2) throw std::runtime_error("cannot reopen input FASTQ files");
            while (read_fastq(r1, rec1)) {
                if (!read_fastq(r2, rec2)) throw std::runtime_error("mate FASTQ ended during output pass");
                if (duplicate[pair_index++]) continue;
                const size_t len1 = bwa_trim_length(rec1.quality, opt.trim_quality);
                const size_t len2 = bwa_trim_length(rec2.quality, opt.trim_quality);
                const bool keep1 = len1 >= static_cast<size_t>(opt.min_length);
                const bool keep2 = len2 >= static_cast<size_t>(opt.min_length);
                if (keep1 && keep2) {
                    write_fastq(out1, rec1, len1);
                    write_fastq(out2, rec2, len2);
                    ++kept_pairs;
                    quality_bases += len1 + len2;
                } else if (keep1) {
                    write_fastq(singleton, rec1, len1);
                    ++singleton_reads;
                    quality_bases += len1;
                } else if (keep2) {
                    write_fastq(singleton, rec2, len2);
                    ++singleton_reads;
                    quality_bases += len2;
                } else {
                    ++short_pairs;
                }
            }
        }

        std::ofstream report(opt.report);
        if (!report) throw std::runtime_error("cannot open report file");
        report << "metric\tvalue\n"
               << "total_input_pairs\t" << total_pairs << '\n'
               << "dedup_quality_eligible_pairs\t" << quality_eligible_pairs << '\n'
               << "duplicate_pairs_removed\t" << duplicate_pairs.load() << '\n'
               << "paired_reads_after_trimming\t" << kept_pairs * 2 << '\n'
               << "kept_pairs_after_trimming\t" << kept_pairs << '\n'
               << "singleton_reads_after_trimming\t" << singleton_reads << '\n'
               << "pairs_with_both_mates_short\t" << short_pairs << '\n'
               << "quality_bases_after_trimming\t" << quality_bases << '\n'
               << "min_identical_bases\t" << opt.seed_length << '\n'
               << "max_diff_rate\t" << opt.max_diff_rate << '\n'
               << "min_mean_quality_for_dedup\t" << opt.min_mean_quality << '\n'
               << "bwa_trim_quality_threshold\t" << opt.trim_quality << '\n'
               << "minimum_trimmed_length\t" << opt.min_length << '\n';
        return 0;
    } catch (const std::exception &error) {
        std::cerr << "ERROR: " << error.what() << '\n';
        return 2;
    }
}
