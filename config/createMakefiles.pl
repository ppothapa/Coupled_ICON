#! /usr/bin/env perl
#
# Usage: createMakefiles.pl
#

use Cwd;
use File::Copy;
use Getopt::Long;
use File::Path;

# Option processing

GetOptions( 
	    'target=s'  => \$target,
            'srcdirs=s' => \$srcdirs,
            'with_ocean=s' => \$with_ocean,
	    ) or die "\n\nUsage: config/createMakefiles.pl --target=<OS_CPU> --srcdirs=< list of src directories>\n";  



$target =~ s/\s+//g;

$srcdirs =~ s/[\"\']*//g;
$srcdirs =~ s/\s+/ /g;
$srcdirs =~ s/^\s+//;

(@directories) = split / /, $srcdirs;

# determine base directory

$prefix = &cwd;

# make architecture dependend build directories

print "\n";
print "createMakefiles:\n\n";

$build_path = &BuildSetup ($prefix, $target, \@directories);

# collect include files in build tree

copy ("config/config.h", "${build_path}/include/config.h");

opendir(DIR, "include") or die "Unable to open include:$!\n";
@incs = grep /\.(inc|h)/, readdir(DIR);
closedir(DIR);
foreach $inc ( @incs ) {
    copy ( "include/${inc}", "${build_path}/include/${inc}" );
}

if ( -d "externals/mtime/include" ) {
    opendir(DIR, "externals/mtime/include");
    @incs = grep /\.(inc|h)/, readdir(DIR);
    closedir(DIR);
    foreach $inc ( @incs ) {
	copy ( "externals/mtime/include/${inc}", "${build_path}/include/${inc}" );
    }
}

# scan dependencies (recursive)

foreach $dir ( @directories ) {

# global variables

    %vpath_directories  = ();
    
    %target_programs    = ();
    
    %module_definitions = ();
    %module_usage       = ();
    
    %fortran_includes   = ();
    %c_includes         = (); 
    
    %ifdefs             = ();
    
    @source_files       = ();

    &ScanDirectory ($dir, $dir, 0);

    $print_path = $build_path;
    $print_path =~ s/$prefix\///;

    print "creating $print_path/$dir/Makefile\n";

    mkpath("$build_path/$dir");
    open MAKEFILE, ">$build_path/$dir/Makefile";

    print MAKEFILE "#----------------------------------------------------------\n";
    print MAKEFILE "# Makefile automatically generated by createMakefiles.pl\n";
    print MAKEFILE "#----------------------------------------------------------\n";
    print MAKEFILE "\n";

    $add_vpath_level = 0;
    if ( "$dir" ne "src" ) {
	if ( $dir =~ m/^externals/) {
	    @subdirs = split(/\//, $dir);
	    print MAKEFILE "SHELL = /bin/bash\n\n";
	    print MAKEFILE "LIB  = $subdirs[1]\n\n";
	    $add_vpath_level = 2;
	} else {
	    print MAKEFILE "LIB  = $dir\n\n";
	}
    }

# write VPATH
    
    @vpath = ();
    push @vpath, "VPATH = ";
    while ( my ($key, $value) = each(%vpath_directories) ) {
	if ( $dir ne "src" ) { $value++; }
#	for my $i ( 0 .. $value ) {
#	    #$key = "../".$key;
#            # RS Hack to work with JSBACH sub-directories in src/lnd_phy_jsbach
#            $key = "../".$key unless $key =~ /^..\/..\/..\/src\/lnd_phy_jsbach\/.*/ ;
#	}
        # Use a constant upward path, this allows arbitary source folder tree depth 
        $key = "../../../".$key.":";
	if ($add_vpath_level == 2) {
	    $key = "../../".$key;
	}
	push @vpath, $key;
    }
    print MAKEFILE @vpath;
    print MAKEFILE "\n\n";
    
# write compile and link information

    print MAKEFILE "%.o: %.f\n";
    print MAKEFILE "\t\$(F77) \$(F77FLAGS) -c \$<\n";
    print MAKEFILE "\n";

    if (($dir =~ m/^externals/) && ($target =~ /^sx/)) {
	print MAKEFILE "\n";
	print MAKEFILE "%.o: %.f90\n";
	print MAKEFILE "\t\$(FC) \$(FlibFLAGS) -c \$<\n";
	print MAKEFILE "\n\n";
    } else {	
	print MAKEFILE "%.o: %.f90\n";
	print MAKEFILE "\t\$(FC) \$(FFLAGS) -c \$<\n";
	print MAKEFILE "\n";
    }

#     print MAKEFILE "%.obj: %.f90\n";
#     print MAKEFILE "\t\$(FC) \$(FFLAGS) -c \$<\n";
#     print MAKEFILE "\n";
  

# write all source files but not the program files
    
    %seen = ();
    @seen{values(%target_programs)} = ();
    
    @sources = ();    
    foreach my $file (@source_files) {
	push (@sources, $file) unless exists $seen{$file};
    }

    if ( "$dir" eq "support" ) {
	($cpu, $vendor, $os) = split /-/, $target;
	if ( "$os" eq "superux") {
	    push (@sources, "rtc_sx.s");
	}
    } 

    print MAKEFILE "SRCS =\t";
    &PrintWords(8, 0, @sources);
    print MAKEFILE "\n\n";
    
    print MAKEFILE "OBJS =\t";
    @objects = ();
    foreach $src (@sources) {
	$src =~ s/\.\S+$/\.o/;
	push @objects, $src;
    }
    &PrintWords(8, 0, @objects);
    print MAKEFILE "\n\n";
    
# write targets

    @target_all = ();
    while ( my ($key, $value) = each(%target_programs) ) {
	push @target_all, $key;
    }
    if ( "$dir" ne "src" ) {
	if ( $dir =~ m/^externals/) {
	    print MAKEFILE ".PHONY: \$(LIB)\n\n";
	    print MAKEFILE "all: \$(LIB)\n\n";
	    print MAKEFILE "\$(LIB): ../../../lib/lib\$(LIB).a\n";
	} else {
	    print MAKEFILE "all: \$(LIB)\n\n";
	    print MAKEFILE "\$(LIB): ../lib/lib\$(LIB).a\n";
	}
    } else {
	print MAKEFILE "all: \$(OBJS) ";
	&PrintWords (13, 0, @target_all);
    }
    print MAKEFILE "\n\n";

    if ( "$dir" ne "src" ) {
	if ( $dir =~ m/^externals/) {
	    print MAKEFILE "../../../lib/lib\$(LIB).a: \$(OBJS)\n";
	    print MAKEFILE "\t\$(AR) \$(ARFLAGS) ../../../lib/lib\$(LIB).a \$(OBJS)\n";
            print MAKEFILE "\t\@for modfile in \$(wildcard *.mod); do \\\n";
            print MAKEFILE "\t\tcp \$\$modfile ../../../include; \\\n"; 
            print MAKEFILE "\t done\n\n";
	    $include_dir = $dir;
	    $include_dir =~ s/src/include/;
            print MAKEFILE "CFLAGS += -I../../../../../$include_dir\n";
            print MAKEFILE "FFLAGS := \$(subst ../module,../../../module, \$(FFLAGS))\n";	    
            print MAKEFILE "\n\n";
	} else {
	    print MAKEFILE "../lib/lib\$(LIB).a: \$(OBJS)\n";
	    print MAKEFILE "\t\$(AR) \$(ARFLAGS) ../lib/lib\$(LIB).a \$(OBJS)\n\n";
	}

	if ( "$dir" eq "support" ) {
            print MAKEFILE "ifeq (\$(ARCH), SX)\n";
            print MAKEFILE "rtc_sx.o: rtc_sx.s\n";
	    print MAKEFILE "\t\$(AS) -c rtc_sx.s\n";
	    print MAKEFILE "endif\n\n";
	}

    } else {
	print MAKEFILE ".PHONY: version.c\n\n";
	print MAKEFILE "version.c:\n";
	print MAKEFILE "\t../../../config/pvcs.pl --srcdir ../../..\n\n";
	print MAKEFILE "version.o: version.c\n\n";
	while ( my ($key, $value) = each(%target_programs) ) {
	    my $okey = $key;
	    $okey =~ s/ *$/.o/;	
	    print MAKEFILE "$okey: $value\n";
	    print MAKEFILE "$key: $okey \$(OBJS) version.o\n";
	    print MAKEFILE "\t\$(FC) \$(LDFLAGS) -o ../bin/\$@ \$< \$(OBJS) version.o \$(LIBS)\n";
	}
    }
    
    print MAKEFILE "clean:\n";
    print MAKEFILE "\trm -f *.o ../module/*.mod\n";
    print MAKEFILE "\n";
    
# print Fortran module dependencies

# don't need c implicit rules are used and don't like dependecies on system header files
#    @c_sources = grep /\.c$/, @sources;
#
#    if ( $#sources > 0 ) {
#        for my $file ( @c_sources ) {
#	    my ($object) = $file;
#	    $object =~ s/c$/o/;
#	    my (@includes) = @{$c_includes{$file}};
#	    print MAKEFILE "$object: $file ";
#	    &PrintWords (length($object)+length($file)+3, 0, @includes);
#	    print MAKEFILE "\n\n";
#	}
#    }
    
    for my $file ( keys %module_usage ) {
	next if $file =~ /c$/;
	my ($object) = $file;
	$object =~ s/f90$/o/;
	my (@modules) = @{$module_usage{$file}};
	my (@dependencies) = ();
	for $i ( 0 .. $#modules) {		 
	    my ($ofile) = $module_definitions{$modules[$i]};
	    $ofile =~ s/f90$/o/;
	    next if $object =~ $ofile;
	    push @dependencies, $ofile;
	}
	next if $object =~ $file;
	print MAKEFILE "$object: $file ";
	&PrintWords (length($object)+length($file)+3, 0, @dependencies);
	print MAKEFILE "\n\n";
    }
    
    close (MAKEFILE);
    
}
print "\n";
#-----------------------------------------------------------------------------
exit;
#-----------------------------------------------------------------------------

#=============================================================================
#
# Subroutines:
#
#=============================================================================
sub ScanDirectory {

# arguments to function    
    my ($workdir) = shift; 
    my ($workpath) = shift;
    my ($level)   = shift;    

# local work
    my($startdir) = &cwd;

    $level++;

    chdir($workdir) or die "Unable to enter dir $workdir:$!\n";

    my($pwd) = &cwd;
    $pwd =~ s/$prefix//;
    $pwd =~ s/^\///;

#     $vpath_directories{$pwd} = $level; it does not work with symbolic links
    $vpath_directories{$workpath} = $level;
    opendir(DIR, ".") or die "Unable to open $workdir:$!\n";
    my @names = readdir(DIR);
    closedir(DIR);
    
    foreach my $name (@names){
        next if ($name eq "."); 
        next if ($name eq "..");
        next if ($name eq ".svn");
        next if ($name eq "templates");
        next if ($name eq "hydro");
        next if ($name eq "interface");
        next if ($name eq "nh");
        next if ($name eq "phys");
        next if ($name eq "sw_options");
        next if (($with_ocean eq "no") and (($name eq "ocean") or ($name eq "sea_ice")) );

        if (-d $name){
	    $nextpath="$workpath/$name";
            &ScanDirectory($name, $nextpath, $level);
            next;
        } else {
	    if ($name =~ /\.[c|f|F]{1}(90|95|03)?$/) {
		push @source_files, $name;

		open F, "<$name";
		while (<F>) {
		    if (/^ *MODULE/i && ! /procedure/i) {
			s/MODULE//i;
			s/\s//g;
			(my @d0) = split '!';
			$module_definitions{$d0[0]} = $name;
		    }
		    if ( (/^ *PROGRAM/i) && ($name !~ /\.c$/) ) {
			s/PROGRAM//i;
			s/\s//g;
			(my @d1) = split '!';
			$target_programs{$d1[0]} = $name;
		    }
		    if (/^ *USE/i) {
			s/USE//i;
			(my @d2) = split ',';
			$d2[0] =~ s/\s//g;
			if (exists $module_usage{$name}) {
			    push @{ $module_usage{$name} }, $d2[0]; 
			} else {
			    $module_usage{$name}[0] = $d2[0];
			}
		    }
		    if (/^ *INCLUDE/i) {
			(my @d3) = split '\'';
			$d3[1] =~ s/\s//g;
			if (exists $fortran_includes{$name}) {
			    push @{ $fortran_includes{$name} }, $d3[1]; 
			} else {
			    $fortran_includes{$name}[0] = $d3[1];
			}
		    }
		    if (/^ *#include/) {
			if ( $name =~ /\.c$/ ) {
			    (my @d4) = split ' ';
			    $d4[1] =~ s/\s//g;
			    if (exists $c_includes{$name}) {
				push @{ $c_includes{$name} }, $d4[1]; 
			    } else {
				$c_includes{$name}[0] = $d4[1];
			    }
			} elsif ( $name =~ /\.[f|F]{1}(90|95|03)?$/ ) {
			    (my @d4) = split ' ';
			    $d4[1] =~ s/\s//g;
			    if (exists $fortran_includes{$name}) {
				push @{ $fortran_includes{$name} }, $d4[1]; 
			    } else {
				$fortran_includes{$name}[0] = $d4[1];
			    }
			}
		    }
		}
		close (F);
	    }
	}
    }
    chdir($startdir) or die "Unable to change to dir $startdir:$!\n";
}
#-----------------------------------------------------------------------------
sub PrintWords {

# arguments to function
    my ($columns) = 78 - shift(@_);
    my ($extratab) = shift(@_);
    my ($wordlength);
    
# local work
    print MAKEFILE @_[0];
    $columns -= length(shift(@_));
    foreach my $word (@_) {
	$wordlength = length($word);
	if ($wordlength + 1 < $columns) {
	    print MAKEFILE " $word";
	    $columns -= $wordlength + 1;
	}
	else {
	    if ($extratab) {
		print MAKEFILE " \\\n\t\t$word";
		$columns = 62 - $wordlength;
            }
	    else {
		print MAKEFILE " \\\n\t$word";
		$columns = 70 - $wordlength;
            }
	}
    }
}
#-----------------------------------------------------------------------------
sub BuildSetup {
    
# arguments to function
    my ($build_path) = shift;
    my ($build_target) = shift;
    my ($build_directories) = shift;

# local work
    my @path7 = ();

    my ($path1) = "$build_path/build";
    my ($path2) = "$path1/$build_target";
    my ($path3) = "$path2/include";
    my ($path4) = "";
    if ( $build_target !~ /superux/ ) {
	$path4 = "$path2/module";
    }
    my ($path5) = "$path2/bin";
    my ($path6) = "$path2/lib";
    foreach $dir ( @{$build_directories} ) {    
	push @path7, $path2."/".$dir;
    }

    if ( -d $build_path ) {
	if ( ! -d $path1 ) {
	    mkdir $path1 || die "Couldn't create build directory";
	}
	mkdir "$path2" || die "Couldn't create build directory";
    } else {
        die "Couldn't create build directory";
    }
    
    if ( ! -d $path3 ) {
	mkdir $path3 || die "Couldn't create module directory";		
    }
    if ( (! -d $path4) && ($path4 =~ /module/) ) {
	mkdir $path4 || die "Couldn't create bin directory";		
    }
    if ( ! -d $path5 ) {
	mkdir $path5 || die "Couldn't create lib directory";		
    }
    if ( ! -d $path6 ) {
	mkdir $path6 || die "Couldn't create lib directory";		
    }
    foreach $path ( @path7 ) {        
	if ( ! -d $path ) {
	    mkdir $path || die "Couldn't create object directory";
	}		
    }

    return $path2;
}
#=============================================================================


