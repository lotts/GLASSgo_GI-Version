#!/usr/bin/perl
# version 1.2
#
# MIT License
#
# Copyright (c) 2017 
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:

# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.

# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.

use Getopt::Std;
use Env;
use File::Basename;
use File::Spec;
use strict;
use Newick;
use DistanceFactory;

use Bio::Matrix::IO;
use Bio::AlignIO;
use Bio::TreeIO;


# analyze user input and execute the corresponding workflow
my %opts;
getopts('f:n:d:o:t:m:h', \%opts);
my ($operation, $output_path, $threshold, $method) = checkUserInput(%opts);

# MAIN OPERATION - WORKFLOWS
if ($operation eq "help") {
    # call help function
    help();
    exit;
}elsif ($operation eq "exit"){
    # exit program
    exit;
}elsif($operation eq "f"){
    fasta2lumberjack($opts{"f"}, $output_path, $threshold, $method);
}elsif($operation eq "n"){
    newick2lumberjack($opts{"n"}, $output_path, $threshold, $method);
}elsif($operation eq "d"){
    distMatrix2lumberjack($opts{"d"}, $output_path, $threshold, $method);
}

# =====================================================================


# ---------------------------------------------------------------------
# FUNCTIONS - check user input
# ---------------------------------------------------------------------

=head2 checkDependencies

 Title   : checkDependencies
 Usage   : my $boolean = checkDependencies();
 Function: The function analyses the path variables and looking for RNApdist.
           If it is defined in the path variables, the function returns true otherwise false.
 Returns : "true" or "false"
 Args    : -

=cut

sub checkDependencies{
    my $tool_name = "RNApdist";
    my $tool_path = '';
    for my $path ( split /:/, $ENV{PATH} ) {
        if ( -f "$path/$tool_name" && -x _ ) {
            $tool_path = "$path/$tool_name";
            last;
        }
    }
    if ($tool_path eq '') {
        print("RNApdist cannot be found. Please install RNApdist\n\t(http://www.tbi.univie.ac.at/RNA/RNApdist.html)\nand make it available as environmental variable!\n");
        return "false";
    }else{
        return "true";
    }
}

=head2 checkUserInput

 Title   : checkUserInput
 Usage   : my $args = checkUserInput(%opts);
 Function: Evaluates the user input and decides which algorithm will be used or help function
           will be called. =
 Returns : (str, str/int, float) - (algorithm_type, file_path, distance_threshold)
 Args    : %opts includes command line input -> getopts('f:n:d:o:t:h', \%opts);

=cut

sub checkUserInput{
    my %opts = @_;
    my $algorithm, my $threshold, my $method;
    # -f fasta_folder      -n newick string    -t threshold (2 default)
    # -d distance matrix   -o out path         -m split method (1 default) -h help file
    #getopts('f:n:d:o:t:h', \%opts);
    if (exists $opts{"h"}) {
        return "help",0,0;
    }
    my $num_para = 0;
    my $count = 0;
    foreach my $i (sort keys %opts){
        if ($i ne "h" && $i ne "o" && $i ne "t" && $i ne "m") {
            $num_para += 1;
            $algorithm = $i;
        }
    }
    # ERROR - some dependencies are missing e.g. RNApdist
    if ($algorithm eq "f" && checkDependencies() eq "false") {
        return "exit",0,0;
    }
    # ERROR - to many arguments
    if ($num_para == 0) {
        print("To less arguments, please use either '-f or -n or -d' to specify the input!\n\t=> londen -h to get more informations\n");
        return "exit",0,0;
    }
    # ERROR - to many arguments
    if ($num_para > 2) {
        print("To many arguments, please use either '-f or -n or -d' to specify the input!\n\t=> londen -h to get more informations\n");
        return "exit",0,0;
    }
    # check parameter "output path"
    my $out_path = "";
    if (exists $opts{"o"}) {
        $out_path = $opts{"o"};
    }else{
        #system("mkdir londen_result");
        $out_path = "./londen_result/";
    }
    # check parameter "threshold"
    if (exists $opts{"t"}) {
        if ($opts{"t"} > 0) {
            $threshold = $opts{"t"};
        }else{
            $threshold = 2;
        }
    }else{
        $threshold = 2;
    }
    
    # check parameter "method" 1:=single cut; 2:=double cut
    if (exists $opts{"m"}) {
        if (($opts{"m"} == 1) or ($opts{"m"} == 2)) {
            $method = $opts{"m"};
        }else{
            $method = 1;
        }
    }else{
        $method = 1;
    }
    
    return $algorithm, $out_path, $threshold, $method;
}

=head2 help

 Title   : help
 Usage   : help();
 Function: Prints the usage of Lumberjack.
 Returns : -
 Args    : -

=cut

sub help{
    print("usage londen:" . "\n");
    print("arguments:\n\t-f fasta file or path to fasta file(s)\n\t-n newick file or path to newick file(s)\n\t-d distance file or path to distance matrix file(s)\n\t-o path for results (default: ./londen_result/)\n");
    print("\t-t threshold to detect inconsistent edges (default: 2)\n");
    print("example-1:\n\tlonden -f ./myFastaFiles/ -o ./myLondenResults/\n");
    print("example-2:\n\tlonden -f ./myFastaFiles/my.fasta\n");
}


# ---------------------------------------------------------------------
# WORKFLOW - compute clusters from fasta files (complete pipeline)
# ---------------------------------------------------------------------

=head2 fasta2lumberjack

 Title   : fasta2lumberjack
 Usage   : fasta2lumberjack(path_to_fasta, out_path, distance_threshold);
 Function: The workflow works with fasta file(s) and consists of 7 steps.
           Step 1: execute RNApdist from vienna package
           Step 2: compute neighbor joining tree
           Step 3: compute adjacency matrix, internal nodes and external nodes
           Step 4: analyze tree and remove inconsistent nodes
           Step 5: prepare new clusters for output
           Step 6: store clusters
           Step 7: remove tmp files
 Returns : -
 Args    : str, str, float (>0)

=cut

sub fasta2lumberjack{
    my($in_path, $out_path, $threshold, $method) = @_;
    #my $newick_out = "outputree_tmp.nw";
    my $newick_out = $$;
    my @seqFiles = ();
     
    if (-d $in_path) {
        opendir(IMD, $in_path) || die("Cannot open fasta directory!");
        @seqFiles = readdir(IMD);
        closedir(IMD);
    }elsif(-e $in_path){
        push(@seqFiles, $in_path);
    }else{
        help();
        exit;
    }
    
    # iterate through all fasta files
    foreach my $f (@seqFiles){
        unless ( ($f eq ".") || ($f eq "..") ){
            my $fastaSeqs = "", my $switch ="";
            if (-d $in_path) {
                my $abs_path   = File::Spec->rel2abs($in_path);
                $fastaSeqs  = $abs_path . "/" . $f;
                $switch = "true";
            }else{
                $fastaSeqs  = $f;
                $switch = "false";
            }
            
            # execute RNApdist from vienna package
            my $resultRNApdist = computeRNApdist($fastaSeqs);
            my ($dist_path) = rnaPdist2matrix($resultRNApdist);            
            
            # compute neighbor joining tree
            my $dfactory = DistanceFactory->new(-method => 'NJ');
            my $tmp_newick_out = ">" . $newick_out;
            my $treeout = Bio::TreeIO->new(-format => 'newick', -file => $tmp_newick_out);
            my $parser = Bio::Matrix::IO->new(-format => 'phylip', -file => $dist_path);
            my $mat  = $parser->next_matrix;
            my $tree = $dfactory->make_tree($mat);
            $treeout->write_tree($tree);
            
            # compute adjacency matrix, internal nodes and external nodes
            my ($adjacencyMatrix, $internalNodes, $externalNodes) = computeAdjacencyPlusNodes($newick_out);
            
            # analyze tree and remove inconsistent nodes
            findAndRemoveInconsistentEdges($internalNodes, $adjacencyMatrix, $threshold, $method);
            
            # prepare new clusters for output
            my $clusterOut = adjaceny2clusters($adjacencyMatrix, $externalNodes);
            
            # store clusters 
            saveCluster2File($f, $fastaSeqs, $clusterOut, $out_path, $switch);
            
            # remove temp files
            system("rm " . $newick_out);
            system("rm " . $dist_path);
        }
    }
}

# ---------------------------------------------------------------------
# WORKFLOW - compute clusters from newick file
# ---------------------------------------------------------------------

=head2 newick2lumberjack

 Title   : newick2lumberjack
 Usage   : newick2lumberjack(path_to_fasta, out_path, distance_threshold);
 Function: The workflow works with newick file(s) and consists of 4 steps.
           Step 1: compute adjacency matrix, internal nodes and external nodes
           Step 2: analyze tree and remove inconsistent nodes
           Step 3: prepare new clusters for output
           Step 4: store clusters
 Returns : -
 Args    : str, str, float (>0)

=cut

sub newick2lumberjack{
    my($in_path, $out_path, $threshold, $method) = @_;
    my @seqFiles = ();
     
    if (-d $in_path) {
        opendir(IMD, $in_path) || die("Cannot open fasta directory!");
        @seqFiles = readdir(IMD);
        closedir(IMD);
    }elsif(-e $in_path){
        push(@seqFiles, $in_path);
    }else{
        help();
        exit;
    }
    
    # iterate through all newick files
    foreach my $f (@seqFiles){
        unless ( ($f eq ".") || ($f eq "..") ){
            my $newickSeqs = "", my $switch ="";
            if (-d $in_path) {
                my $abs_path = File::Spec->rel2abs($in_path);
                $newickSeqs  = $abs_path . "/" . $f;
                $switch = "true";
            }else{
                $newickSeqs  = $f;
                $switch = "false";
            }
            
            # compute adjacency matrix, internal nodes and external nodes
            my ($adjacencyMatrix, $internalNodes, $externalNodes) = computeAdjacencyPlusNodes($newickSeqs);
            
            # analyze tree and remove inconsistent nodes
            findAndRemoveInconsistentEdges($internalNodes, $adjacencyMatrix, $threshold, $method);
            
            # prepare new clusters for output
            my $clusterOut = adjaceny2clusters($adjacencyMatrix, $externalNodes);
            
            # store clusters 
            saveNewickCluster2File($f, $clusterOut, $out_path, $newickSeqs);
        }
    }
}



# ---------------------------------------------------------------------
# WORKFLOW - compute clusters from distance-matrix file
# ---------------------------------------------------------------------

=head2 distMatrix2lumberjack

 Title   : distMatrix2lumberjack
 Usage   : distMatrix2lumberjack(path_to_fasta, out_path, distance_threshold);
 Function: The workflow works with fasta file(s) and consists of 7 steps.
           Step 1: compute neighbor joining tree
           Step 2: compute adjacency matrix, internal nodes and external nodes
           Step 3: analyze tree and remove inconsistent nodes
           Step 4: prepare new clusters for output
           Step 5: store clusters
           Step 6: remove tmp files
 Returns : -
 Args    : str, str, float (>0)

=cut

sub distMatrix2lumberjack{
    my($in_path, $out_path, $threshold, $method) = @_;
    #my $newick_out = "outputree_tmp.nw";
    my $newick_out = $$;
    my @seqFiles = ();
    
    if (-d $in_path) {
        opendir(IMD, $in_path) || die("Cannot open fasta directory!");
        @seqFiles = readdir(IMD);
        closedir(IMD);
    }elsif(-e $in_path){
        push(@seqFiles, $in_path);
    }else{
        help();
        exit;
    }
    
    # iterate through all distance matrices  files
    foreach my $f (@seqFiles){
        unless ( ($f eq ".") || ($f eq "..") ){ 
            my $dist_path = "", my $switch ="";
            if (-d $in_path) {
                my $abs_path   = File::Spec->rel2abs($in_path);
                $dist_path  = $abs_path . "/" . $f;
                $switch = "true";
            }else{
                $dist_path  = $f;
                $switch = "false";
            }
            
            # compute neigbhor joining tree
            #my $dfactory = Bio::Tree::DistanceFactory->new(-method => 'NJ');
            my $dfactory = DistanceFactory->new(-method => 'NJ');
            my $tmp_newick_out = ">" . $newick_out;
            my $treeout = Bio::TreeIO->new(-format => 'newick', -file => $tmp_newick_out);
            my $parser = Bio::Matrix::IO->new(-format => 'phylip', -file => $dist_path);
            my $mat  = $parser->next_matrix;
            my $tree = $dfactory->make_tree($mat);
            $treeout->write_tree($tree);
            
            # compute adjacency matrix, internal nodes and external nodes
            my ($adjacencyMatrix, $internalNodes, $externalNodes) = computeAdjacencyPlusNodes($newick_out);
            
            # analyze tree and remove inconsistent nodes
            findAndRemoveInconsistentEdges($internalNodes, $adjacencyMatrix, $threshold, $method);
            
            # prepare new clusters for output
            my $clusterOut = adjaceny2clusters($adjacencyMatrix, $externalNodes);
            
            # store clusters 
            saveNewickCluster2File($f, $clusterOut, $out_path, $newick_out, $switch);
            
            # remove temp files
            system("rm " . $newick_out);
        }
    }
}





# ---------------------------------------------------------------------
# ---------------------------------------------------------------------
# FUNCTIONS
# ---------------------------------------------------------------------

=head2 computeRNApdist

 Title   : computeRNApdist
 Usage   : my $rna_pdist_dist_matrix = computeRNApdist($fasta_file);
 Function: Executes RNApdist with a specific fasta file and returns a distance matrix.
 Returns : <distance_str>RNApdist
 Args    : str

=cut

sub computeRNApdist{
	my $fastaSeqs = $_[0];
	## execute RNApdist from vienna package - make the tool available in your path variables!
	my $resultRNApdist = `RNApdist -Xm < $fastaSeqs`;
	## clean output files because they are not needed for the further steps
	system("rm .\/*.ps");
    #print ($resultRNApdist . "\n");
	return $resultRNApdist;
}

=head2 rnaPdist2matrix

 Title   : rnaPdist2matrix
 Usage   : my $phylip_file = rnaPdist2matrix($rna_pdist_dist_matrix);
 Function: Transform distance string to matrix (array).
 Returns : path to phylip matrix
 Args    : <distance_str>RNApdist

=cut

sub rnaPdist2matrix{
	my $resultRNApdist = $_[0];
	my @oneDArr        = (); # includes fasta identifier
	my @matrixLines    = ();
    my $filename = "./foo_dist_matrix.tmp";
	# convert "$resultRNApdist" into a matrix
	my @tmpRes = split('\n', $resultRNApdist);
	
	foreach (@tmpRes){	
		if($_ =~ m/^>\S/){
			push(@oneDArr, $_);
		}
		if($_ =~ m/^\d/ || $_ =~ m/^-\d/){
            #print($_ . "\n");
			push(@matrixLines, $_);
		}
	}
	my @arrayMatrix = ();
	#init diagonal from '@arrayMatrix' with 'x'
	for(my $i = 0 ; $i < scalar(@oneDArr) ; $i++){
		$arrayMatrix[$i][$i] = 0.0;
	}
	foreach (@matrixLines){	
		my @tmp = split(' ', $_);
		for(my $i = 0 ; $i < scalar(@tmp) ; $i++){
			$arrayMatrix[$i][scalar(@tmp)] = $tmp[$i];
			$arrayMatrix[scalar(@tmp)][$i] = $tmp[$i];
		}
	}
    # open tmp_file
    open(my $fh, '>', $filename) or die "Could not open file '$filename' $!";
    print $fh "\t" . scalar(@oneDArr) . "\n";
    # iterate through $arrayMatrix
    for(my $i = 0 ; $i < scalar(@oneDArr) ; $i++){
        print $fh $oneDArr[$i] . "\t";
        for(my $j = 0 ; $j < scalar(@oneDArr) ; $j++){
            print $fh $arrayMatrix[$j][$i] . " ";
        }
        print $fh "\n";
    }
    close $fh;
    # create phylip matrix and store on hdd (return file name)
	return $filename;
}

=head2 computeAdjacencyPlusNodes

 Title   : computeAdjacencyPlusNodes
 Usage   : my ($adjacencyMatrix, $internalNodes, $externalNodes) = computeAdjacencyPlusNodes();
 Function: Reads a file in Newick format and generates an adjacency matrix,
           hash of inner nodes and external nodes (leafs).
 Returns : \@adjacencyMatrix, \%internalNodes, \%externalNodes;
 Args    : $path_newick - path to newick file

=cut

sub computeAdjacencyPlusNodes{
	my $tmpFilename =  $_[0];
	my @adjacencyMatrix = (); 
	# read from file and load the newick stored file from above
	my $tree = Newick->from_file($tmpFilename);
    
	my $maxRange = 0;
	# init adjacency matrix with zeros
	for my $node (1..$tree->nodes-1) {
		my $parent = $tree->parent->[$node];
		if($maxRange < $node){
			$maxRange = $node;
		}
		if($maxRange < $parent){
			$maxRange = $parent;
		}
	}
	for my $i (0..$maxRange) {
		for my $j (0..$maxRange) {
			$adjacencyMatrix[$i][$j] = -1;
		}	
	}
	my %externalNodes = (); # external nodes are leafs
	my %internalNodes = (); # internal nodes could be connected to leafs or to another internal nodes
	
	for my $node (1..$tree->nodes-1) {
		my $parent = $tree->parent->[$node];
		my $tmpEdge = $node . "-" . $parent;
		my $tmpDist = abs($tree->branch_length->[$node]);
		# collect all connections and nodes
		if( $tree->node_name->[$node] ne "" ){
			# add leafs
			$externalNodes{$node}{"parentNode"} = $parent;
			$externalNodes{$node}{"nodeName"}   = $tree->node_name->[$node];
			$internalNodes{$parent}{"externalEdges"}{$tmpEdge} = $tmpDist;
		}else{
			# collect internal nodes + all direct neighborhood connections
			$internalNodes{$node}{"internalEdges"}{$tmpEdge}   = $tmpDist; 
			$internalNodes{$parent}{"internalEdges"}{$tmpEdge} = $tmpDist;
		}
		# compute adjacency matrix
		$adjacencyMatrix[$node][$parent] = $tmpDist;
		$adjacencyMatrix[$parent][$node] = $tmpDist;
	}
	return \@adjacencyMatrix, \%internalNodes, \%externalNodes;
}

=head2 findAndRemoveInconsistentEdges

 Title   : findAndRemoveInconsistentEdges
 Usage   : findAndRemoveInconsistentEdges($internalNodes, $adjacencyMatrix, $threshold);
 Function: The function analyses the adjacency matrix, taking into account local threshold
           and internal nodes list. The algorithm computes local distances and removes
           inconsistence edges from the adjacency matrix
 Returns : works on references
 Args    : $internalNodes (%hash), $adjacencyMatrix (@array), $threshold (float)

=cut

sub findAndRemoveInconsistentEdges{
	my $hashRef = $_[0];
	my $arrRef  = $_[1]; # adjacency matrix
    my $threshold = $_[2];
    my $method = $_[3];
	my @updateList = ();

    
    my %tmp_hash = ();
	# compute sum of all outgoing edges
	for my $key ( keys %$hashRef ) {
		my $sum     = 0;
		my $members = 0;
		
		for my $keyA ( keys %{$hashRef->{$key}->{"internalEdges"}} ) {
			if(exists $hashRef->{$key}->{"internalEdges"}->{$keyA}){
				$sum += $hashRef->{$key}->{"internalEdges"}->{$keyA};
				$members++;
			}
		}
		for my $keyA ( keys %{$hashRef->{$key}->{"externalEdges"}} ) {
			if(exists $hashRef->{$key}->{"externalEdges"}->{$keyA}){
				$sum += $hashRef->{$key}->{"externalEdges"}->{$keyA};
				$members++;
			}
		}
        # ******* double cut ********
        if ($method == 2) {
            # mode 2 - NEW - check every outgoing internal edge, if l > 2*mean -> delete connection
            for my $keyA ( keys %{$hashRef->{$key}->{"internalEdges"}} ) {
                if(exists $hashRef->{$key}->{"internalEdges"}->{$keyA}){
                    # compute threshold
                    my $tresh = $threshold * ($sum / $members);
                    my @tmpEdge = split('\-', $keyA);
                    # update adjacency matrix
                    if($hashRef->{$key}->{"internalEdges"}->{$keyA} > $tresh){
                        if(exists $tmp_hash{$tmpEdge[0]}{$tmpEdge[1]}){
                            $tmp_hash{$tmpEdge[0]}{$tmpEdge[1]} += -1;
                        }else{
                            $tmp_hash{$tmpEdge[0]}{$tmpEdge[1]} = -1;
                        }
                    }
                }
            }
        }
        
		# check every outgoing internal edge, if l > 2*mean -> delete connection
		for my $keyA ( keys %{$hashRef->{$key}->{"internalEdges"}} ) {
			if(exists $hashRef->{$key}->{"internalEdges"}->{$keyA}){
				# compute threshold
				my $tresh = $threshold * ($sum / $members);
				if($hashRef->{$key}->{"internalEdges"}->{$keyA} > $tresh){
					my @tmpEdge = split('\-', $keyA);
					# update adjacency matrix
                    if ($method == 1) {
                        my $pairOne = $tmpEdge[1] . "-" . $tmpEdge[0];
                        push(@updateList, $keyA);
                        push(@updateList, $pairOne);
                        $arrRef->[$tmpEdge[0]]->[$tmpEdge[1]] = -1;
                        $arrRef->[$tmpEdge[1]]->[$tmpEdge[0]] = -1;
                    }elsif(($tmp_hash{$tmpEdge[0]}{$tmpEdge[1]}) < -1){
                        my $pairOne = $tmpEdge[1] . "-" . $tmpEdge[0];
                        push(@updateList, $keyA);
                        push(@updateList, $pairOne);
                        $arrRef->[$tmpEdge[0]]->[$tmpEdge[1]] = -1;
                        $arrRef->[$tmpEdge[1]]->[$tmpEdge[0]] = -1;
                    }
				}
			}
		}
	}
	# update all nodes (connections to delete)
	for my $key ( keys %$hashRef ) {
		for my $edge (0..@updateList-1){
			if(exists $hashRef->{$key}->{"internalEdges"}->{$updateList[$edge]}){
				delete $hashRef->{$key}->{"internalEdges"}->{$updateList[$edge]};
			}	
		}
	}
}

=head2 saveCluster2File

 Title   : saveCluster2File
 Usage   : saveCluster2File($file_name, $fasta_path, $clusters ,$out_path, $toggle);
 Function: Creates new fasta files corresponding to the new clusters and stores clusters on hdd.
 Returns : -
 Args    : $file_name (string), $fasta_path (string), $clusters (%hash), $out_path (string),
           $toggle "true", "false" ("true" for directory as input, "false" for single file as input)

=cut

sub saveCluster2File{
    my $filename   = $_[0];
	my $fastaPath  = $_[1];
	my $clusterOut = $_[2];
	my $outPath    = $_[3];
    my $switch = $_[4]; # true = iterate of directory, false = only one file 
    
	# get file name 
	my($fastaIdent, $pathToSave, $suffix) = fileparse($fastaPath);
	my @clusterIdentifier = split('-', $fastaIdent);
	
	# read cluster fasta file
	my %currFastaHash = ();
    my $startSeq = 0;
    my $newSeq = 0;
    my $tmpKey;
	open(FILE, "<$fastaPath"); #$fastaSeqs only for testing , change variable into $seqs
		while(<FILE>){
			chomp($_);
			if($startSeq && !($_ eq "") && !($_ =~ m/^>/)){
				$currFastaHash{$tmpKey} .= $_;	
			}
			if($newSeq && !($_ eq "")){
				$currFastaHash{$tmpKey} = $_;  		
		  		$newSeq   = 0;
		  		$startSeq = 1;
			}
			if($_ =~ m/^>/){
				$tmpKey = $_;
				$currFastaHash{$tmpKey} = "";
				$newSeq   = 1;
				$startSeq = 0;
			}
		}
	close(FILE);
	
	my %valiKeys = ();
	for my $key ( keys %$clusterOut ) {
		for my $keyB ( keys %{$clusterOut->{$key}} ) {
			if(exists $currFastaHash{ $clusterOut->{$key}->{$keyB} }){
				$valiKeys{$key} = 1;
				last;
			}
		}
	}
		
	my $clusterSize = keys( %valiKeys );
	my $currNum     = 0;
	for my $key ( keys %valiKeys ) {
		$currNum++;
		
		my $tmpOut = "";
		$tmpOut = "cluster-" . $currNum . "of" . $clusterSize . "-";
		
		my @tmpSplit = split('-', $clusterIdentifier[-1]);
		my $diffReads   = 0;
		my %diffGenomes = ();
		for my $keyB ( keys %{$clusterOut->{$key}} ) {
			if(exists $currFastaHash{ $clusterOut->{$key}->{$keyB} }){
				my @tmpGenome = split(/-|\+/, $clusterOut->{$key}->{$keyB});
				if(exists $diffGenomes{$tmpGenome[0]}){
					$diffGenomes{$tmpGenome[0]} +=1
				}else{
					$diffGenomes{$tmpGenome[0]}  =1
				}
				$diffReads++;
			}
		}
        
        my $tmpOutput = "";
        if ($switch eq "true") {
            $tmpOutput = $outPath . $tmpOut . $filename;
        }else{
            my @tmpName = split("\/", $filename);
            $tmpOutput = $outPath . $tmpOut . $tmpName[-1];
        }
        
		print "CLUSTER\n";
		for my $keyB ( keys %{$clusterOut->{$key}} ) {
            if(exists $currFastaHash{ $clusterOut->{$key}->{$keyB} }){
                print $clusterOut->{$key}->{$keyB} . "\n";
                print $currFastaHash{$clusterOut->{$key}->{$keyB}} . "\n";
            }
        }
	}	
}

=head2 saveNewickCluster2File

 Title   : saveNewickCluster2File
 Usage   : saveNewickCluster2File($file_name_newick, $clusters ,$out_path, $newick_path, $toggle);
 Function: Creates new cluster files corresponding to the new "clusters" ($clusters) and stores clusters on hdd.
 Returns : -
 Args    : $file_name_newick (string), $clusters (%hash), $out_path (string), $newick_path (string),
           $toggle "true", "false" ("true" for directory as input, "false" for single file as input)
 

=cut

sub saveNewickCluster2File{
    my $newick_file = $_[0];
	my $clusterOut = $_[1];
	my $outPath    = $_[2];
    my $newick_path = $_[3];
    my $switch = $_[4]; # true = iterate of directory, false = only one file 
    
    my $tree = Newick->from_file($newick_path);
    my %node_names = ();
	# init adjacency matrix with zeros
	for my $node (1..$tree->nodes-1) {
        if (defined $tree->node_name->[$node]) {
            $node_names{$tree->node_name->[$node]} = 1;
        }
    }
    my $tmpOut = "";
    if ($switch eq "true") {
        $tmpOut = $outPath . "londen_" . $newick_file;
    }else{
        my @tmpName = split("\/", $newick_file);
        $tmpOut = $outPath . "londen_" . $tmpName[-1];
    }

	my %outClusters = ();
	for my $key ( keys %$clusterOut ) {
		for my $keyB ( keys %{$clusterOut->{$key}} ) {  
            if (exists $node_names{$clusterOut->{$key}->{$keyB}}) {
                if (exists $outClusters{$key}) {
                    $outClusters{$key} .= ";" .$clusterOut->{$key}->{$keyB};
                }else{
                    $outClusters{$key} = $clusterOut->{$key}->{$keyB};
                }
            }
		}
	}
        
    open(FILE, ">$tmpOut");
    my $cl = 1;
    foreach my $key ( keys %outClusters ) {
        print FILE "cluster " . $cl . ": " . $outClusters{$key} . "\n";
        $cl += 1;
    }
	close(FILE);
}

=head2 adjaceny2clusters

 Title   : adjaceny2clusters
 Usage   : my $clusters = adjaceny2clusters($adjacencyMatrix, $externalNodes);
 Function: Transfers adjacency matrix into a list of clusters. This function will be
           applied after findAndRemoveInconsistentEdges() was executed.
 Returns : \%clusters
 Args    : $adjacencyMatrix (@array), $externalNodes (%hash)

=cut

sub adjaceny2clusters{
	my $arrRef        = $_[0]; 	# adjacency matrix
	my $externalNodes = $_[1];  # -> 
								# $externalNodes{$node}{"parentNode"} = $parent;
								# $externalNodes{$node}{"nodeName"}   = $tree->node_name->[$node];
	my %nodes = ();		my %clusters = ();		my @heap = ();
	
	# (1) generate hash with keys as number of matrix rows
	for my $i (0..@$arrRef-1){
		$nodes{$i} = 1;	
	}
	# (2) start at node 0 and follow the path
	# if path is finished, but there exists another nodes in
    # the nodes-hash then pick some nodes from the hash and compute another path
	# an path corresponds to one cluster
	my $curr      = 0;
	my %wasInHeap = ();
	
	while ( keys( %nodes ) > 0 ) {
		# every while loop corresponds to one cluster , one graph
		my $tk =  (keys %nodes)[-1];
		push(@heap, $tk);
		
		while((scalar(@heap) > 0)){
			$tk = pop(@heap);
			if(exists $externalNodes->{$tk}->{"nodeName"}){
				$clusters{$curr}{$tk} = $externalNodes->{$tk}->{"nodeName"};
			}else{
				$clusters{$curr}{$tk} = 1;
			}
			delete $nodes{$tk};			
			
			for my $i (0..@$arrRef-1){
				if($arrRef->[$tk]->[$i] != -1){
					if(exists $wasInHeap{$i}){
					}else{
						$wasInHeap{$i} = 1;
						push(@heap, $i);
					}
				}
			}
		}
		$curr++;
	}
	# return %cluster as a reference, 
	# %clusters{clusterNumber}{nodeNumber}=>if "leaf"" than "leafName" else "1"
	return \%clusters;
}
