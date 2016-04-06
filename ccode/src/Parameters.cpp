#include "Parameters.hpp"
#include <iostream>

Parameters::Parameters() : graph_filename(""), c(0), rank(0),
                           num_iter(0), oversampling(0), U_filename(""), D_filename("") {}

bool Parameters::init(int argc, char* argv[]) {
    if (argc < 15) {
        printf("Current number of arguments = %d and is insufficient\n", argc);
        exit_with_help();
        return false;
    }
    for (int i = 1; i < argc; i++) {
        if (argv[i][0] != '-') {
            std::cout << "The keys have to start with -" << std::endl;
            return false;
        }
        if (++i >= argc) {
            exit_with_help();
            return false;
        }
        switch (argv[i-1][1]) {
            case 'g':
                graph_filename = argv[i];
                break;

            case 'r':
                rank = atoi(argv[i]);
                break;

            case 'c':
                c = atof(argv[i]);
                break;

            case 'i':
                num_iter = atoi(argv[i]);
                break;
                
            case 's':
                oversampling = atoi(argv[i]);
                break;

            case 'u':
                U_filename = argv[i];
                break;
            
            case 'd':
                D_filename = argv[i];
                break;
                
            default:
                fprintf(stderr, "Unknown option: -%c\n", argv[i-1][1]);
                exit_with_help();
                return false;
        }
    }
    return true;
}

void Parameters::exit_with_help() {
    printf(
    "Usage: simrank_lowrank [options]\n"
    "options:\n"
    "-g graph file name\n"
    "-r rank of approximation\n"
    "-c SimRank parameter\n"
    "-u filename where to save the lowrank matrix U\n"
    "-d filename where to save the diagonal of the diagonal matrix D\n"
    "-i number of iteration in lowrank approximation\n"
    "-s oversampling parameter"
    );
}
