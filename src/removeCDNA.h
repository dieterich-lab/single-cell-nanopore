#include <seqan/sequence.h>
#include <seqan/arg_parse.h>
#include <seqan/seq_io.h>

#include <iostream>
#include <seqan/align.h>
#include <seqan/bam_io.h>
#include<iostream>
#include<fstream>
#include <seqan/find.h>
#include <stdlib.h>     /* srand, rand */
#include <time.h>       /* time */
#include <map>

using namespace seqan;
using namespace std;

/*

void writeRead(SeqFileOut & seqFileOut, string & last_id, Dna5String & end1, Dna5String & end2, bool rc, bool verbose = false)
{
    if(rc){
        Dna5StringReverseComplement rcend1(end1);
        Dna5StringReverseComplement rcend2(end2);

        if(verbose){
            std::cout << "Use reverse Complement of reads to restore original orientation.\n";
            std::cout << "before: " << end1 << "\n";
            std::cout << "before2: " << end2 << "\n";
        }
        Dna5String tmpend1 = rcend1;
        end1 = rcend2;
        end2 = tmpend1;
        if(verbose){
            std::cout << "after:  " << end1 << "\n";
            std::cout << "after2:  " << end2 << "\n";
        }
    }
    else
    {
        if(verbose)
            std::cout << "Is in original orientation\n";
    }

    CharString idend1 = last_id;
    CharString idend2 = last_id;
    idend1 += "_end1";
    idend2 += "_end2";

    writeRecord(seqFileOut, idend1, end1);
    writeRecord(seqFileOut, idend2, end2);
}
*/

bool checkIfSameOrientation(Dna5String & originalRead, Dna5String & bamRead){

    bool reverseComplement;
    Dna5StringReverseComplement rcBamRead(bamRead);
//     appendValue(rcReads, myModifier);
    for(int i = 0; i < length(bamRead); ++i)
    {
//         std::cout << i << ": " << originalRead[i] << "\t" << bamRead[i] << "\t" << rcBamRead[i] << "\n";

        if(originalRead[i] != bamRead[i] && originalRead[i] != rcBamRead[i]){
            std::cout << "Bam and fasta reads are not the same.\n";
            exit(0);
        }
        if(originalRead[i] != bamRead[i]){
            reverseComplement = true;
            break;
        }
        if(originalRead[i] != rcBamRead[i]){
            reverseComplement = false;
            break;
        }
    }
/*
    if(reverseComplement)
        std::cout << "reverseComplement\n";
    else
        std::cout << "not\n";*/

    return reverseComplement;
}

Dna5String & getOriginalRead(auto & readMap, string & readID)
{
    CharString creadID = readID;
    auto search = readMap.find(creadID);
    if (search == readMap.end())
    {
        std::cout << "Did not find " << readID << " in fasta file." << "\n";
        exit(0);
    }
    Dna5String & originalRead = search->second;

    return originalRead;
}

void removeCDNA(Options &o, std::vector<BamAlignmentRecord > & records)
{
    /*
    ArgumentParser parser("Search");

    addOption(parser, ArgParseOption("b", "bam", "Path to the bam file", ArgParseArgument::INPUT_FILE, "IN"));
    setRequired(parser, "bam");

    addOption(parser, ArgParseOption("f", "fasta", "Path to the fasta with original reads to recover orientation", ArgParseArgument::INPUT_FILE, "IN"));
    hideOption(getOption(parser, "fasta"));

    addOption(parser, ArgParseOption("o", "output", "Path to fasta output file", ArgParseArgument::OUTPUT_FILE, "OUT"));
    setRequired(parser, "output");

    addOption(parser, ArgParseOption("p", "first", "First p reads", ArgParseArgument::INTEGER, "INT"));
    hideOption(getOption(parser, "first"));

    addOption(parser, ArgParseOption("mt", "threshold", "Number of matches required to define start or end of the cDNA part", ArgParseArgument::INTEGER, "INT"));

    addOption(parser, ArgParseOption("v", "verbose", ""));

    ArgumentParser::ParseResult res = parse(parser, argc, argv);
    if (res != ArgumentParser::PARSE_OK)
        return res == ArgumentParser::PARSE_ERROR;

    CharString bamPath, outputPath;
    CharString fastaPath = "";
    int batchSize1 = 100000;

    int barcode_umi_length = 30;
    int threshold = 5;
    int first = 9999999;

    getOptionValue(bamPath, parser, "bam");
    getOptionValue(fastaPath, parser, "fasta");
    getOptionValue(outputPath, parser, "output");
//     getOptionValue(barcodeLength, parser, "barcodeL");
    getOptionValue(first, parser, "first");
    getOptionValue(threshold, parser, "threshold");
    bool verbose = isSet(parser, "verbose");*/
/*
    //Load original reads

    //prepare Bam
    BamFileIn bamFileIn;
    BamHeader header;
    if (!open(bamFileIn, toCString(bamPath)))
    {
        std::cerr << "ERROR: could not open input file " << bamPath << ".\n";
        return 1;
    }
    readHeader(header, bamFileIn);

    SeqFileOut seqFileOut(toCString(outputPath));*/

    //TODO add threshold as parameter
    int threshold = 5;

    //modify bam (Add UMI and Barcode (determined from flexbar) to each Alignment record)
    CigarElement<> soft = seqan::CigarElement<>('S', 1);
    CigarElement<> match = seqan::CigarElement<>('M', 1);
    CigarElement<> del = seqan::CigarElement<>('D', 1);
    auto op_soft = soft.operation;
    auto op_match = match.operation;
    auto op_del = del.operation;

    bool verbose = false;

//     while (!atEnd(bamFileIn))
    for(uint32_t i = 0; i < records.size(); ++i)
    {
        auto & record = records[i];
            Dna5String end1 = record.seq;
            Dna5String end2 = record.seq;
            CharString id = record.qName;
            CharString qual1 = record.qual;
            CharString qual2 = record.qual;


            if(verbose){
                std::cout << "Accession in Bam " << toCString(id) << "\n";
                std::cout << "Read length: " << length(record.seq) << ": " << record.seq << "\n";
            }

            String<CigarElement<> > cigar = record.cigar;
            //take front
            int pos = 0;
            //count amount of overall matches in cigar string
            int match_count = 0;
            bool found = false;
            for(int i = 0; i < length(cigar); ++i){

                if (cigar[i].operation == op_match)
                    match_count += cigar[i].count;

                //taake prefix of the read as soon as more than threshold many matches are found
                if (match_count > threshold && !found){
//                     std::cout << "found pos: " << pos << "\n";
                    end1 = prefix(end1, pos + 1);
                    qual1 = prefix(qual1, pos + 1);
                    found = true;
                }
                //calculate end position of prefix according to cigar string
                if(!found && (cigar[i].operation == op_match || cigar[i].operation == op_soft || cigar[i].operation == op_del))
                    pos += cigar[i].count;
            }

            int score = 0;


            if(!found){
                if(verbose)
                    std::cout << "Take whole read since it was not aligned\n";
                o.fastaRecords.push_back(make_tuple(id, end1, qual1));
//                 std::make_tuple(10, "Test", 3.14, std::ref(n), n);
                continue;
            }

            pos = 0;
            int match_count2 = 0;
            bool found2 = false;
            for(int i = length(cigar) - 1; i >= 0 ; --i){

                if (cigar[i].operation == op_match)
                    match_count2 += cigar[i].count;

                //since the overall match count is know we can stop as soon as we determined the start position of the prefix
                if (match_count2 > threshold){
    //                     std::cout << "found pos: " << pos << "\n";
                    end2 = suffix(end2, length(end2) - pos - 1);
                    qual2 = suffix(qual2, length(qual2) - pos - 1);
                    found2 = true;
                    break;
                }

                if(cigar[i].operation == op_match || cigar[i].operation == op_soft || cigar[i].operation == op_del)
                    pos += cigar[i].count;
            }

                if(verbose){
                    std::cout << "lenght end1: " << length(end1) << "\n";

                    std::cout << "lenght end2: " << length(end2) << "\n";

                    std::cout << "Score: " << (match_count) << "\n";
                }


                score = match_count;


                bool rc = hasFlagRC(record);

                if(rc){
                    Dna5StringReverseComplement rcend1(end1);
                    Dna5StringReverseComplement rcend2(end2);
                    ModifiedString<CharString, ModReverse> rcqual1(qual1);
                    ModifiedString<CharString, ModReverse> rcqual2(qual2);

                    if(verbose){
                        std::cout << "Use reverse Complement of reads to restore original orientation.\n";
                        std::cout << "before: " << end1 << "\n";
                        std::cout << "before2: " << end2 << "\n";
                    }
                    Dna5String tmpend1 = rcend1;
                    CharString tmpqual1 = rcqual1;

                    qual1 = rcqual2;
                    qual2 = tmpqual1;

                    end1 = rcend2;
                    end2 = tmpend1;
                    if(verbose){
                        std::cout << "after:  " << end1 << "\n";
                        std::cout << "after2:  " << end2 << "\n";
                    }
                }
                else
                {
                    if(verbose)
                        std::cout << "Is in original orientation\n";
                }

                CharString idend1 = id;
                CharString idend2 = id;
                idend1 += "_end1";
                idend2 += "_end2";

//                 std::cout << end1 << "\n" << qual1 << "\n";
                assert(length(end1) == length(qual1));
                assert(length(end2) == length(qual2));

                o.fastaRecords.push_back(make_tuple(idend1, end1, qual1));
                o.fastaRecords.push_back(make_tuple(idend2, end2, qual2));

    }


    std::cout << "Finished removing cDNA!\n";
}
