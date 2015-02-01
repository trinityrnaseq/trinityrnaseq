#include <iostream>
#include <string>
#include <fstream>
#include <cstdlib>
#include <sstream>
#include <omp.h>
#include <vector>
#include <map>
#include <algorithm>
#include <signal.h>
#include <sys/wait.h>

#include "argProcessor.hpp"


string usage();


using namespace std;


static int VERBOSE_LEVEL = 0;

int main(int argc, char* argv[]) {
    
    /* This program executes  commands in Parallel using C++, the program is a wrapper designed for the orignal Trinity Program
       for seemless parallel OPENMP execution, Authors of this wrapper are MB Couger (mbcouger(AT Symbol)gmail.com, Matt Stowe mstowe(AT Symbol)okstate.edu*/ 
    // code modified by bhaas for general cmd-line processing


    ArgProcessor args(argc, argv);
    
    bool shuffle_flag = false;


    if (args.isArgSet("-help") 
        || ! (args.isArgSet("-c") ) 
        || ! (args.isArgSet("-CPU") ) 
        ) {
        cerr << usage() << endl << endl;
        return(1);
    }
    
    if (args.isArgSet("-v")) {
        VERBOSE_LEVEL = 1;
    }
    else if (args.isArgSet("-vv")) {
        VERBOSE_LEVEL = 2;
    }

    if (args.isArgSet("-shuffle")) {
        shuffle_flag = true;
    }
    

    string commands_file = args.getStringVal("-c");
    
    int num_cpus = args.getIntVal("-CPU");
        
    
    string failed_commands_filename = "FailedCommands";
    if (args.isArgSet("-failed_cmds")) {
        failed_commands_filename = args.getStringVal("-failed_cmds");
    }
    
    //declare variables start input stream
    ifstream in;    // Create an input file stream.
    in.open(commands_file.c_str());  // Use it to read from a file named data.txt.
    if ( ! in ) {
        stringstream errstr;
        errstr << "Error, cannot open file : " << commands_file << endl;
        cerr << errstr.str();
        exit(1);
    }
    

    // check to see if any commands had completed
    map<string,bool> previously_completed_commands;
    bool have_previously_completed_commands = false;
    
    ofstream successfully_completed_fh;

    {
        ifstream prev_completed_fh;
        stringstream completed_commands_filename_ss;
        completed_commands_filename_ss << commands_file << ".completed";
        string completed_commands_filename = completed_commands_filename_ss.str();
        prev_completed_fh.open(completed_commands_filename.c_str());
        
        if (prev_completed_fh) {
            have_previously_completed_commands = true;
            
            string command_line;
            getline(prev_completed_fh, command_line);
            while(! prev_completed_fh.eof()) {
                previously_completed_commands[command_line] = true;
                getline(prev_completed_fh, command_line);
            }
            prev_completed_fh.close();
        }

        successfully_completed_fh.open(completed_commands_filename.c_str(), ios_base::app); // open for append
    }
    
    long int NumberofCommands;
    NumberofCommands=0;
    
    vector<string> commandsArray;

    string line;
    getline(in,line);
    
    while (! in.eof()) {

        if (have_previously_completed_commands 
            &&
            previously_completed_commands.find(line) != previously_completed_commands.end()) {
            
            // this command has been run before successfully, not running it again.
            if (VERBOSE_LEVEL) {
                cerr << "warning, command: " << line << " has successfully completed from a previous run.  Skipping it here." << endl;
            }
            
            
        }
        else {
            NumberofCommands++;
            commandsArray.push_back(line);
        }
        getline(in,line);
    }
    
    cerr << "Number of  Commands: " <<  NumberofCommands << endl;
           
    //Parrell Execution of Individual Commands 
    vector<string> failedCommands;
    int num_failed_commands = 0;
    int num_succeeded_commands = 0;

    if (shuffle_flag) {
        random_shuffle(commandsArray.begin(), commandsArray.end());
    }

    omp_set_dynamic(0);
    omp_set_num_threads(num_cpus);

    #pragma omp parallel for schedule (dynamic)
    for (long int i=0;i<NumberofCommands;i++) {
        string command = commandsArray[i];         
        
        if (VERBOSE_LEVEL == 2) {

            int thread_no = omp_get_thread_num();

            #pragma omp critical (standard_error)
            {
                cerr << "CMD[" << i << "], thread[" << thread_no << "]: " << command << endl;
            }
        }
        
        int ret = system(command.c_str());

        // exit if child received SIGINT or SIGQUIT
        if (WIFSIGNALED(ret) &&
            (WTERMSIG(ret) == SIGINT || WTERMSIG(ret) == SIGQUIT)) {
            #pragma omp critical (exit_critical)
            exit(ret); // behavior undefined if exit() called more than once
        }
        else if (ret != 0) {
                        
            #pragma omp critical (capture_failed_command)
            {
                num_failed_commands++;
                failedCommands.push_back(command);
            }
            

            if (VERBOSE_LEVEL == 2) {
                #pragma omp critical (standard_error)
                {
                    cerr << "FAILURE:[" << i << "]  " << command << endl;
                }
            }

        }
        else {
            #pragma omp critical (report_success)
            {
                num_succeeded_commands++;
            
                successfully_completed_fh << command << endl;
            }
                
            if (VERBOSE_LEVEL == 2) {
            #pragma omp critical (standard_error)
                cerr << "SUCCESS:[" << i << "]  " << command << endl;
            }
        }
        
        if (VERBOSE_LEVEL == 1) {
            
            stringstream ss;
            ss << "\rsucceeded(" << num_succeeded_commands << ")";
            if (num_failed_commands > 0) {
                ss << ", failed(" << num_failed_commands << ")";
            }
            int total_executed = num_succeeded_commands + num_failed_commands;
            float percent_done = (float)total_executed/NumberofCommands * 100;
            
            ss << "   " << percent_done << "% completed.    ";

            #pragma omp critical (standard_error)
            cerr << ss.str();
        }
    }
    
    successfully_completed_fh.close();
    

    //ErrorCheckingOut
    if (num_failed_commands != 0) {
        ofstream outdata(failed_commands_filename.c_str());
        for (int t=0;t<num_failed_commands; ++t) {
            outdata << failedCommands[t]  << endl;
        }
        outdata.close();
        cout << endl << endl << "We are sorry, commands in file: [" << failed_commands_filename << "] failed.  :-( " << endl << endl;
        
        exit(1);
    }
    else {
        cout << endl << endl << "All commands completed successfully. :-)" << endl << endl;
        exit(0); // used to be return(0), but sometimes in OMP land this would not exit 0....?!?!
    }
    
}

string usage () {

    stringstream ss;

    ss << endl 
       << "##########################################################" << endl
       << "#" << endl
       << "# Usage: ParaFly (opts)" << endl
       << "#" << endl
       << "# Required: " << endl
       << "#   -c <str>              :filename containing list of bash-style commands to execute." << endl
       << "#   -CPU <int>            :number_of_threads" << endl
       << "#" << endl
       << "# Optional:" << endl
       << "#   -shuffle              :randomly shuffles the command order. " << endl
       << "#   -failed_cmds <str>    :filename to capture failed commands.  default(\"FailedCommands\")" << endl
       << "#   -v                    :simple progress monitoring." << endl
       << "#   -vv                   :increased verbosity in progress monitoring." << endl
       << "#" << endl
       << "##########################################################" << endl << endl;
    
    ss << "Note: This process creates a file named based on your commands filename with a .completed extension." << endl
       << "This enables a resume functionality, where if rerun, only those commands not completed successfully will be reprocessed." << endl << endl;

    return(ss.str());
}

