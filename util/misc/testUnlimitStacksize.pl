#!/usr/bin/env perl

use strict;
use warnings;
no strict 'subs';
#use BSD::Resource;

main: {
    
    print "\n\nBefore unlimit stack space:";
    system("bash -c ulimit -a");
    print "Attempting to unlimit stack space...";
    try_unlimit();

    

    
    print "After unlimit stack space: ";
    system("bash -c ulimit -a");
    system("ls");
    my $ls = `ls`;
    print "$ls\n";
    
    exit(0);
}


sub try_unlimit {
    #eval "
     #           use BSD::Resource;
      #  setrlimit(RLIMIT_STACK, RLIM_INFINITY, RLIM_INFINITY);";


    my $unset_stacksize_code = "use BSD::Resource;"
        . "setrlimit(RLIMIT_STACK, RLIM_INFINITY, RLIM_INFINITY);"
        . "my (\$soft, \$hard) = getrlimit(RLIMIT_STACK);"
        . "print \"stack_soft: \$soft, stack_hard: \$hard  \";"
        . "if (\$soft != -1) { die \"Couldn't unset stacksize\";}";


    eval($unset_stacksize_code);
    
    if( $@ ) {
        warn <<"EOF";
        
        $@
                
        Unable to set unlimited stack size. Please install the BSD::Resource
            Perl module to allow this script to set the stack size, or set it
            yourself in your shell before running Trinity (ignore this warning if
            you have set the stack limit in your shell). See the following URL for
            more information:

            http://trinityrnaseq.sourceforge.net/trinity_faq.html#ques_E

EOF
;
    }
    else {
    
            
        print "Successfully set unlimited stack size.\n";
        print "###################################\n\n";



    }
    
    ## verify
    print "\n\n";
    eval("use BSD::Resource;"
         ."my (\$soft, \$hard) = getrlimit(RLIMIT_STACK);"
         ."print \"Verifying: stack_soft: \$soft, stack_hard: \$hard  \";"  );
    
    

        
    return;
}

