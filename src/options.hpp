#pragma once
#include "annotations.hpp"
#include "store.hpp"
#include <getopt.h>

namespace gc {
struct Options {
    StoreOptions store;
    AnnotationOptions annotation;
    std::string chr_map, missing_contig = "warn";
    bool overwrite = false;
};
inline std::set<std::string> parse_list(const std::string& value, const std::string& option) {
    std::set<std::string> result;
    for (const auto& item : split(value,',')) {
        require(!item.empty() && trim(item) == item, "Invalid " + option + " list: " + value);
        result.insert(item);
    }
    return result;
}
inline Options parse_options(int argc, char** argv, const char* help) {
    enum { EMIT=1000, MAX_DISTANCE, WINDOW, DOT_BP, MAX_DOTS, LABEL, CODING, BIOTYPES,
           EXCLUDE_BIOTYPES, FEATURES, FORMAT, CHR_MAP, MISSING_CONTIG, PREFIX, OVERWRITE,
           SORT_MEM, TEMP_DIR, QUIET, TSS_SOURCE, CANONICAL_TAG, PROMOTER_UP, PROMOTER_DOWN,
           FEATURE_TRANSCRIPTS, TRANSCRIPT_FEATURES };
    static const option opts[] = {
        {"annotation",required_argument,nullptr,'a'}, {"emit",required_argument,nullptr,EMIT},
        {"max-distance",required_argument,nullptr,MAX_DISTANCE}, {"window",required_argument,nullptr,WINDOW},
        {"dot-bp",required_argument,nullptr,DOT_BP}, {"max-dots",required_argument,nullptr,MAX_DOTS},
        {"label",required_argument,nullptr,LABEL}, {"protein-coding",no_argument,nullptr,CODING},
        {"biotypes",required_argument,nullptr,BIOTYPES}, {"exclude-biotypes",required_argument,nullptr,EXCLUDE_BIOTYPES},
        {"features",required_argument,nullptr,FEATURES}, {"annotation-format",required_argument,nullptr,FORMAT},
        {"chr-map",required_argument,nullptr,CHR_MAP}, {"missing-contig",required_argument,nullptr,MISSING_CONTIG},
        {"prefix",required_argument,nullptr,PREFIX}, {"overwrite",no_argument,nullptr,OVERWRITE},
        {"sort-mem",required_argument,nullptr,SORT_MEM}, {"temp-dir",required_argument,nullptr,TEMP_DIR},
        {"tss-source",required_argument,nullptr,TSS_SOURCE}, {"canonical-tag",required_argument,nullptr,CANONICAL_TAG},
        {"promoter-upstream",required_argument,nullptr,PROMOTER_UP}, {"promoter-downstream",required_argument,nullptr,PROMOTER_DOWN},
        {"feature-transcripts",required_argument,nullptr,FEATURE_TRANSCRIPTS},
        {"transcript-features",required_argument,nullptr,TRANSCRIPT_FEATURES},
        {"quiet",no_argument,nullptr,QUIET}, {"help",no_argument,nullptr,'h'}, {nullptr,0,nullptr,0}
    };
    Options result;
    bool emit_given = false, biotypes_given = false;
    optind = 1;
    opterr = 0;
    int c;
    while ((c = getopt_long(argc,argv,"a:h",opts,nullptr)) != -1) {
        const std::string arg = optarg ? optarg : "";
        switch (c) {
        case 'a':
            require(result.store.annotation.empty(), "Supply one --annotation file");
            result.store.annotation = arg;
            break;
        case 'h':
            std::fputs(help,stderr);
            std::exit(EXIT_SUCCESS);
        case EMIT: {
            if (!emit_given) {
                result.annotation.emit.clear();
                emit_given = true;
            }
            auto values = parse_list(arg,"--emit");
            std::set<std::string> allowed;
            for (const auto& field : fields) allowed.insert(field.option);
            for (const auto& value : values) {
                require(value == "all" || allowed.count(value), "Unknown --emit group: " + value);
                if (value == "all") result.annotation.emit.insert(allowed.begin(),allowed.end());
                else result.annotation.emit.insert(value);
            }
            break;
        }
        case MAX_DISTANCE: result.annotation.max_distance = integer(arg,0,INT32_MAX,"--max-distance"); break;
        case WINDOW: result.annotation.window = integer(arg,0,INT32_MAX,"--window"); break;
        case DOT_BP: result.annotation.dot_bp = integer(arg,1,INT32_MAX,"--dot-bp"); break;
        case MAX_DOTS: result.annotation.max_dots = integer(arg,1,1000000,"--max-dots"); break;
        case LABEL:
            require(arg == "name" || arg == "id", "--label must be name or id");
            result.annotation.label_ids = arg == "id";
            break;
        case CODING: case BIOTYPES:
            require(!biotypes_given,"Use only one of --protein-coding and --biotypes");
            biotypes_given = true;
            result.store.filter.include = parse_list(c == CODING ? "protein_coding" : arg,"--biotypes");
            break;
        case EXCLUDE_BIOTYPES: result.store.filter.exclude = parse_list(arg,"--exclude-biotypes"); break;
        case FEATURES: result.store.filter.features = parse_list(arg,"--features"); break;
        case TRANSCRIPT_FEATURES: result.store.filter.transcript_features = parse_list(arg,"--transcript-features"); break;
        case TSS_SOURCE:
            require(arg == "gene" || arg == "canonical","--tss-source must be gene or canonical");
            result.annotation.canonical_tss = arg == "canonical";
            break;
        case CANONICAL_TAG:
            require(!arg.empty() && trim(arg) == arg && arg.find_first_of(",\t\r\n") == arg.npos,
                    "--canonical-tag must be one nonempty transcript tag");
            result.annotation.canonical_tag = arg;
            break;
        case PROMOTER_UP: result.annotation.promoter_upstream = integer(arg,0,INT32_MAX,"--promoter-upstream"); break;
        case PROMOTER_DOWN: result.annotation.promoter_downstream = integer(arg,0,INT32_MAX,"--promoter-downstream"); break;
        case FEATURE_TRANSCRIPTS:
            require(arg == "all" || arg == "canonical","--feature-transcripts must be all or canonical");
            result.annotation.canonical_features = arg == "canonical";
            break;
        case FORMAT:
            require(arg == "auto" || arg == "gtf" || arg == "gff3", "--annotation-format must be auto, gtf or gff3");
            result.store.filter.format = arg == "gtf" ? Format::Gtf : arg == "gff3" ? Format::Gff3 : Format::Auto;
            break;
        case CHR_MAP: result.chr_map = arg; break;
        case MISSING_CONTIG:
            require(arg == "warn" || arg == "error" || arg == "ignore", "--missing-contig must be warn, error or ignore");
            result.missing_contig = arg;
            break;
        case PREFIX:
            require(!arg.empty() && ((arg[0] >= 'A' && arg[0] <= 'Z') || (arg[0] >= 'a' && arg[0] <= 'z') || arg[0] == '_') &&
                    arg.find_first_not_of("abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789_.") == arg.npos,
                    "--prefix must be a VCF tag prefix beginning with a letter or underscore");
            result.annotation.prefix = arg;
            break;
        case OVERWRITE: result.overwrite = true; break;
        case SORT_MEM: result.store.sort_bytes = static_cast<size_t>(integer(arg,1,1048576,"--sort-mem"))*1024*1024; break;
        case TEMP_DIR: result.store.temp_dir = arg; break;
        case QUIET: result.store.quiet = true; break;
        default: throw std::runtime_error("Unknown or incomplete plugin option: " + std::string(argv[optind-1]) + "; use bcftools +genecontext --help");
        }
    }
    require(optind == argc,"Unexpected positional plugin argument; the VCF input must precede --");
    require(!result.store.annotation.empty(),"--annotation is required after --");
    // Feature requirements depend on the final --emit and transcript settings.
    result.store.selection = result.annotation.selection();
    result.store.intron_context = result.annotation.intron_context();
    if (result.store.selection.transcripts) {
        for (const auto& feature : result.store.filter.transcript_features)
            require(!result.store.filter.features.count(feature),"Gene and transcript feature types overlap: "+feature);
    }
    if (result.store.selection.exons || result.store.selection.cds) {
        for (const auto* feature : {"exon","CDS"})
            require(!result.store.filter.features.count(feature) && !result.store.filter.transcript_features.count(feature),
                    std::string("Reserved structural feature type in gene/transcript feature list: ")+feature);
    }
    return result;
}
}
