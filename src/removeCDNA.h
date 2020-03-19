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
        Dna5String end1 = prefix(record.seq, o.fragmentLength);
        Dna5String end2 = suffix(record.seq, length(record.seq) - o.fragmentLength);
        CharString id = record.qName;
        CharString qual1 = prefix(record.qual, o.fragmentLength);
        CharString qual2 = suffix(record.qual, length(record.seq) - o.fragmentLength);


        if(verbose){
            std::cout << "Accession in Bam " << toCString(id) << "\n";
            std::cout << "Read length: " << length(record.seq) << ": " << record.seq << "\n";
        }

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
