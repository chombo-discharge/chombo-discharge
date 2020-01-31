#! /usr/bin/perl

$exe = $ARGV[0];
$procs = $ARGV[1];
if($procs == '') { $procs = 1; }

$totaltime0 = time;

$unamea = `uname -a`;

if($procs <= 16) {
   $nodes = 1;
} elsif ($procs <= 32) {
   $nodes = 2;
} elsif ($procs <= 48) {
   $nodes = 3;
} elsif ($procs <= 64) {
   $nodes = 4;
} else {
   print "Too many procs to run interactively $procs\n";
   exit(0);
}

#$timestamp = `/bin/date '+%m%d%H%M'`;
$date = `date`;

print "running on $procs procs ";
print " ($nodes nodes)   $date\n";
print " system=$unamea \n";

$command = "poe $exe -nodes $nodes -procs $procs  2>eout";
print "$command\n";

$attempt=1;

$unfinished = 1;
while($unfinished) {

   $atime0 = time;
   print "attempt #$attempt : ";
   #system("sleep 1");
   $jobtime0 = time;
   system("$command");
   $jobtime1 = time;
   $commandtime = $jobtime1-$jobtime0;
   
   $unfinished = 0;
   
   open(IN, "eout") || die "can't open eout";
   while (<IN>) {

      #print "$_";
      if (/repo (\S+)/) {
	 $repo = $1;
      } elsif (/Step (\S+) was/) {
	 $stepvalue = $1;
      }

      
      if (/relatively low priority/) {
	 print "interactive queue must be full.  bastards!  ($stepvalue)\n";
	 $unfinished = 1;
	 if( ($attempt+9) % 10 == 0) {
	    print "Step Id          JobName         UserName  Class  ST NDS WallClck Submit Time\n";
	    $llqsi = `llqs | grep 'interac'`;
	    print "$llqsi";
	    #$llqsd = `llqs | grep 'debug'`;
	    #print "$llqsd";

	    #$Njobs = (@jobs = split '\n', $llqs);
	    #system("llqs | grep ' debug '");
	 }
	 break;
      } elsif (/A file or directory in the path name does not exist/) {
	 print "bad filename dumbass!\n";
	 $unfinished = 1;
	 exit(0);
      } elsif (/Starter terminated abnormally/) {
	 print "starter terminated abnormally ??\n";
	 $unfinished = 1;
	 break;
      }
      
   }   
   close(IN);

   $attempt = $attempt+1;
   $atime1 = time;
   $attempttime = $atime1 - $atime0;
   printf("                        %d sec.\n", $attempttime);
}

$user = `whoami`;
print " getnim for user: $user:\n";
system("getnim -U$user");
$totaltime1 = time;
$totaltime = $totaltime1-$totaltime0;
print "done. jobtime: $commandtime sec  total: $totaltime sec,";
print " $attempt attempts (repo used $repo)\n";
