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





void splitReads(Options &o)
{/*
    ArgumentParser parser("Search");

//     addOption(parser, ArgParseOption("b", "bam", "Path to the bam file", ArgParseArgument::INPUT_FILE, "IN"));
//     setRequired(parser, "bam");

    addOption(parser, ArgParseOption("f", "alignedPrimer", "Path to the flexbar P1 primer alignment result fasta", ArgParseArgument::INPUT_FILE, "IN"));
    setRequired(parser, "alignedPrimer");
//     setRequired(parser, "reads1");

    addOption(parser, ArgParseOption("o", "output", "Path to output files prefix", ArgParseArgument::OUTPUT_FILE, "OUT"));
    setRequired(parser, "output");

    addOption(parser, ArgParseOption("c", "checkOrigin", "Check for if Sequence is located on the correct side according to the orientation"));

    addOption(parser, ArgParseOption("l", "readLength", "Shorten longer reads to this length (remove appropiate suffix)", ArgParseArgument::INTEGER, "INT"));

    addOption(parser, ArgParseOption("v", "verbose", ""));

    ArgumentParser::ParseResult res = parse(parser, argc, argv);
    if (res != ArgumentParser::PARSE_OK)
        return res == ArgumentParser::PARSE_ERROR;

    CharString dictPath, flexPath, bamPath, outputPath;
    int wrongSide = 0, forward = 0, reverseC = 0;
    int readLength = 0;

//     getOptionValue(bamPath, parser, "bam");
    getOptionValue(flexPath, parser, "alignedPrimer");
    getOptionValue(outputPath, parser, "output");
    getOptionValue(readLength, parser, "readLength");
//     getOptionValue(batchSize1, parser, "batchSize");
    bool checkOrigin = isSet(parser, "checkOrigin");
    bool verbose = isSet(parser, "verbose");

*/

    //load Primer alignments from fasta


    std::string fastaFile = (o.format == flexbar::FASTA) ? o.targetName + ".fasta" : o.targetName + ".fastq";
    SeqFileIn seqFileInFlex(toCString(fastaFile));


//     CharString outputRevcomp = outputPath;

//     outputPath += "_left_tail_trimmed.fasta";
//     outputRevcomp += "_right_tail_trimmed.fasta";

//     SeqFileOut seqFileOut(toCString(outputPath));
//     SeqFileOut seqFileOutRevcomp(toCString(outputRevcomp));

    int wrongSide = 0, forward = 0, reverseC = 0;
//     int readLength = 0;

    //TODO add as option
    int readLength = 25;
    bool checkOrigin = true;
    bool verbose = false;

    while (!atEnd(seqFileInFlex))
    {
        Dna5String read;
        CharString id;
        CharString qual;
        try
        {

            if(o.format == flexbar::FASTQ)
                readRecord(id, read, qual, seqFileInFlex);
            else
                readRecord(id, read, seqFileInFlex);

            Finder<CharString> finder(id);
            Pattern<CharString, Horspool> pattern("_Flexbar_removal_");
            find(finder, pattern);
            int found = beginPosition(finder);
            //flexbar aligned primer
            if(found > 0){
                Finder<CharString> finder2(id);
                Pattern<CharString, Horspool> pattern("Flexbar_removal_cmdline_rc");
                find(finder2, pattern);
                int end = beginPosition(finder2);
                bool doRC = end > 0;
                bool onRightSide;
                // check If Primer is in correct orientation in regards to aligned cdDNA
                if(checkOrigin){
                    Finder<CharString> finder3(id);
                    Pattern<CharString, Horspool> p_end2("_end2_Flexbar_removal");
                    find(finder3, p_end2);
                    int verify = beginPosition(finder3);
//                     std::cout << toCString(id) << "\n";
                    if (verbose)
                        std::cout << "verify begin Positions:" << verify << "\n";
                    onRightSide = verify > 0;
                }

                if(checkOrigin && ((doRC && !onRightSide) || (!doRC && onRightSide))){
                    ++wrongSide;
                    if(verbose)
                        std::cout << doRC << "\t" << toCString(id) << "\n";
                }
                //flexbar aligned to RC Primer
                else if (doRC){
                    if(readLength > 0 && length(read) > readLength){
                        read = suffix(read, length(read) - readLength);
                        if(o.format == flexbar::FASTQ)
                            qual = suffix(qual, length(qual) - readLength);
                    }
                    ++reverseC;
                    o.rightTail.push_back(make_tuple(id, read, qual));
//                     writeRecord(seqFileOutRevcomp, id, read);
                }
                //flexbar aligned to forward direction of Primer
                else
                {
                    if(readLength > 0 && length(read) > readLength){
                        read = prefix(read, readLength);
                        if(o.format == flexbar::FASTQ)
                            qual = prefix(qual, readLength);
                    }

                    ++forward;
                    o.leftTail.push_back(make_tuple(id, read, qual));
//                     writeRecord(seqFileOut, id, read);
                }
            }
        }
        catch (IOError const & e)
        {
            std::cerr << "ERROR: could not copy record. " << e.what() << "\n";
        }
    }


    close(seqFileInFlex);


    if (checkOrigin){
        std::cout << "Discarded because aligned to wrong side: " << wrongSide << "\n";
        std::cout << "Left: " << forward << "\n";
        std::cout << "Right: " << reverseC << "\n";
    }
    else
    {
        std::cout << "Left: " << forward << "\n";
        std::cout << "Right: " << reverseC << "\n";
    }
//     return 0;
}
