#!/usr/bin/perl -w

use strict;
use warnings;
use Math::Random qw(:all);  # if Perl pukes because modules are not found, check http://www.cpan.org/modules/INSTALL.html
sub trim($);

############################################################################
# mendelian error simulation for genotype vector data on nuclear pedigrees
# written by CVH, cvv8@cornell.edu on or about Sep 19, 2012.
# v2: -simulate phred scores based on exponential distribution
#     -simulate multiple iterations for each parameter configuration
# v3: -change mendelian errors in offspring from 'non-biallelic' model to 'biallelic mendelian incompatible' model
# v4: -add capability to simulate sex linked snps March 1, 2013
# v5: -add mendelian error model to sex chromosomes March 7, 2013
# v6: -change the mendelian error model to 'biallelic' or 'triallelic' variable mendelian error types, also correct missingness to be independent of mendelian error status, i.e. with missingness, some pedigrees could be missing instead of having an error
# to do, a refactor should include subroutines for simulating genotypes(founder,het/homogametic,ParentalGenos) and simError(het/homogametic,sex,ParentalGenos), this moves the logic to the subroutines
############################################################################

#parameters and description
### global variables, do not change these
my @sexLinkage = ();           #sex linked or autosomal transmission
my @errorType = ();            #error model for mendelian error, either biallelic and inconsistent with mendelian transmission, or triallelic where the mendelian error is an offspring homozygous for the third allele
my @MAFS = ();                 #minor allele frequency of the putative biallelic SNP
my @numOffspring = ();         #number of offspring per family
my @numMendErrPerFam = ();     #number of mendelian errors per family, always in the children
my @setGeno = ();              #set the coefficient of variation (as a percent) with this parameter
my @missingness = ();          #zero to one, this proportion of individuals will have no data
my @numFams = ();              #this is the number of pedigrees to simulate for each condition
my $optionLabel = ();          #you guessed it, a label
my $snp = 0;                   #simple counter
#############################################################################


#these are user defined variables, it is ok to change these
my $sexRatio = .5;             #the probability the offspring will be heterogametic.
my $sexUnknown = 0;            #the probability that an offspring with masked (unkonwn) sex
my $heterogameticSexCode = 2;  #code for heterogametic sex, 2 for ZW systems, 1 for XY
my $homogameticSexCode = 1;    #code for homogametic sex, 1 for ZZ systems, 2 for XX, or whatever...
my $genoNull = "";             #null character for genotypes without phred scores, i.e. "" or "."(NOT missing data)
my $genoMissing = "-1";        #symbol for genotypes of individuals who are missing all genotypes
my $iterations = 10000;         #run each of the parameter sets this number of times


open (KEY, ">./simFSJGenotypesv6.key") or die;
print KEY "SNP\tParameter\tParameterValue\tIteration\tParameterSet\n";
open (FAMKEY, ">./simFSJGenotypesv6.famkey") or die;
print FAMKEY "SNP\tFAM\tParameter\tParameterValue\n";
open (OUT, ">./simFSJGenotypesv6.txt") or die;
print OUT "#CHR\tPOS\tSNP\tFAM\tID\tFA\tMO\tSEX\tAA\tAC\tAG\tAT\tCC\tCG\tCT\tGG\tGT\tTT\n";

foreach my $option (8){ #iterate through all of the parameters of interest, or don't by changing this variable
#foreach my $option (1 .. 8){ #iterate through all of the parameters of interest, or don't by changing this variable
#parameter vectors, when the parameter in question is the random variable, this is range of values, else the default value
if ($option == 0){ $optionLabel = "NONE";}
if ($option == 1){ $optionLabel = "MAF"; @MAFS = ("0.01","0.05","0.1","0.25","0.5");} else {@MAFS = (".5");}
if ($option == 2){ $optionLabel = "OFF"; @numOffspring = ("1","2","5","10");} else {@numOffspring = ("2");} #be merciful, keep these parameters low for high iterations
#if ($option == 2){ $optionLabel = "OFF"; @numOffspring = ("1","2","5","10","20","30","40","50","60","70","80","90","100");} else {@numOffspring = ("2");} #be merciful, keep these parameters low for high iterations
if ($option == 3){ $optionLabel = "ERR"; @numMendErrPerFam = ("0","1","2");} else { @numMendErrPerFam = ("0");}
if ($option == 4){ $optionLabel = "GENO";  @setGeno = ("3000","300","1");} else {@setGeno = ("3000");}
if ($option == 5){ $optionLabel = "MIS"; @missingness = ("0","0.1","0.2","0.3","0.4","0.5","0.6","0.7","0.8","0.9","1");} else {@missingness = ("0");}
if ($option == 6){ $optionLabel = "FAM"; @numFams = ("1","10","25","100","500","1000");} else {@numFams = ("20");}
if ($option == 7){ $optionLabel = "SEX"; @sexLinkage = ("auto","sex");} else {@sexLinkage = ("auto");}
if ($option == 8){ $optionLabel = "ETYPE"; @errorType = ("bi","tri");} else {@errorType = ("bi");}

foreach my $maf (@MAFS){
  foreach my $numOff (@numOffspring){
    foreach my $numMendErrPerFam (@numMendErrPerFam){
      foreach my $setGeno (@setGeno){
        foreach my $missingness (@missingness){
          foreach my $numFam (@numFams){
           foreach my $sexLinkage (@sexLinkage){
            foreach my $eType (@errorType){
           foreach my $iteration (1 .. $iterations){
            $snp++;
            #$numFam = 120 / ($numOff + 2); #for set sample sizes
            my $chr = "chr_".$optionLabel."_iter".$iteration."_maf".$maf."_off".$numOff."_fam".$numFam."_err".$numMendErrPerFam."_geno".$setGeno."_mis".$missingness."_sex".$sexLinkage."_eType".$eType;
            print "option:$optionLabel iter:$iteration ParameterSet:$chr\n";
            my $paraValue = '';
            if ($optionLabel eq "NONE")
            {
            	print KEY "$snp\t$optionLabel\tnone\t$iteration\t$chr\n";
            	$paraValue = "none";
            }
            if ($optionLabel eq "MAF")
            {
            	print KEY "$snp\t$optionLabel\t$maf\t$iteration\t$chr\n";
            	$paraValue = $maf;
            }
            elsif ($optionLabel eq "OFF")
            {
            	print KEY "$snp\t$optionLabel\t$numOff\t$iteration\t$chr\n";
            	$paraValue = $numOff;
            }
            elsif ($optionLabel eq "ERR")
            {
            	print KEY "$snp\t$optionLabel\t$numMendErrPerFam\t$iteration\t$chr\n";
            	$paraValue = $numMendErrPerFam;
            }
            elsif ($optionLabel eq "GENO")
            {
            	print KEY "$snp\t$optionLabel\t$setGeno\t$iteration\t$chr\n";
            	$paraValue = $setGeno;
            }
            elsif ($optionLabel eq "MIS")
            {
            	print KEY "$snp\t$optionLabel\t$missingness\t$iteration\t$chr\n";
            	$paraValue = $missingness;
            }
            elsif ($optionLabel eq "SEX")
            {
            	print KEY "$snp\t$optionLabel\t$sexLinkage\t$iteration\t$chr\n";
            	$paraValue = $sexLinkage;
            }
            elsif ($optionLabel eq "ETYPE")
            {
            	print KEY "$snp\t$optionLabel\t$sexLinkage\t$iteration\t$chr\n";
            	$paraValue = $eType;
            }
            else
            {
            	print KEY "$snp\t$optionLabel\t$numFam\t$iteration\t$chr\n";
            	$paraValue = $numFam;
            }

#set up phred score simulation parameterized by exponential random variable
#Aiming for values from actual data!
#highH Mean   : 2833  
#highL Mean   : 326.3  
#lowH Mean   : 302.7  
#lowL Mean   :22.74  
my $lambda1 = $setGeno/10;  #this factor 'default 10' sets the parameter for the exponential distribution for the penultimate most likely genotype as ten times lower (phred score is 10 times lower) than the parameter for the third most likely
my $lambda2 = ($setGeno);
#cludge
my $expOffset = 0; #offset for the phred score distribution, most likely is 0, penultimate is 'expOffset + exp(lambda1)', third is is 'expOffset + exp(lambda2)'
if ($lambda1 == 300){$expOffset = 100;}  #this is cludgy, second offset for high quality genotype configuration

foreach my $fam (1 .. $numFam){  #start the inner most simulation loop for the nuclear family..  bat territory
# print "label $optionLabel chr $chr fam $fam\n";

print FAMKEY "$snp\t$fam\t$optionLabel\t$paraValue\n";

my $famSize = $numOff+2;
my $parent1 = ();
my $parent2 = ();
my $sex = 0;
foreach my $individual (1 .. $famSize){
       my $genotype = ();
       my @genoVector = ($genoNull,$genoNull,$genoNull,$genoNull,$genoNull,$genoNull,$genoNull,$genoNull,$genoNull,$genoNull);
       my $offspringNum = $individual-2; #which offspring, First?  Second? this variable is ignored unless positive
  if ($individual <= 2){ #individual is a parent/founder
       if ($individual == 1){
          if ($sexLinkage eq "sex" && $numMendErrPerFam > 0){  #simulating sex linked snps, force parent1 to be homozygous from Hardy Weinberg proportions
             my $passParentalHet = 0;
             while ($passParentalHet == 0){  #force parent two to be a non-het, but still simulate from the specified population allele frequency, Hardy Weinberg proportions
                  $parent1 = simParent($maf);
                  if ($parent1 ne "12"){
                     $passParentalHet = 1;
                  }
             }
          } else {
          $parent1 = simParent($maf);
          }
        $genotype = $parent1;
       }

       if ($individual == 2){
          my $passParentalHets = 1;
          if ($numMendErrPerFam > 0){ $passParentalHets = 0;}
               $parent2 = simParent($maf);
          if ($parent1 eq "12"){  #if parent one is a het...
             while ($passParentalHets == 0){  #force parent two to be a non-het, but still simulate from the specified population allele frequency, Hardy Weinberg proportions
                  $parent2 = simParent($maf);
                  if ($parent2 ne "12"){
                     $passParentalHets = 1;
                  }
             }
          }
          if ($sexLinkage eq "sex"){  #simulating sex linked snps, force parent2 to be homozygous. Not Hardy Weinberg proportions  In truth, this parent has only one allele of the diplo-genotype.
#          if ($numMendErrPerFam > 0){ die "mendelian errors on sex chromosomes has not been implemented\n";}
            $parent2 = simHeterogameticParent($maf);
          }
          $genotype = $parent2;
       }
       if ($genotype == 11){ $genoVector[0] = 0; $genoVector[1] =sprintf("%.0f",$expOffset+random_exponential(1,$lambda1)); $genoVector[4] = sprintf("%.0f",$expOffset+random_exponential(1,$lambda2));}
       if ($genotype == 12){ $genoVector[1] = 0; $genoVector[0] =sprintf("%.0f",$expOffset+random_exponential(1,$lambda1)); $genoVector[4] = sprintf("%.0f",$expOffset+random_exponential(1,$lambda1));}
       if ($genotype == 22){ $genoVector[4] = 0; $genoVector[1] =sprintf("%.0f",$expOffset+random_exponential(1,$lambda1)); $genoVector[0] = sprintf("%.0f",$expOffset+random_exponential(1,$lambda2));}
       if (rand() < $missingness){ #set the genotypes for this parent to missing
         @genoVector = ($genoMissing,$genoMissing,$genoMissing,$genoMissing,$genoMissing,$genoMissing,$genoMissing,$genoMissing,$genoMissing,$genoMissing);
       }

  }else{ #individual is an OFFSPRING, possibly heterogametic or homogametic, possibly a mendelian error or not, possibly not an error but missing
         #at the moment, sex chromosomes can be missing, but not errors
        if (rand() < $sexRatio){$sex=$heterogameticSexCode;}else{$sex=$homogameticSexCode;} #what about UNKNOWN?
        if ($sexLinkage eq "sex"){
           #this is where to check against a parameterized sex ratio, and determine the sex of the offspring and chromosomal configuration
          if ($sex == $heterogameticSexCode){$genotype = simHeterogameticOffspring($parent1,$parent2);}else{$genotype = simOffspring($parent1,$parent2);}
        } else { #offspring, autosomal snp
          $genotype = simOffspring($parent1,$parent2);
        }
        if ($offspringNum <= $numMendErrPerFam){ #simulate a mendelian error in this offspring
           @genoVector = ($genoNull,$genoNull,$genoNull,$genoNull,$genoNull,$genoNull,$genoNull,$genoNull,$genoNull,$genoNull);
           if ($eType eq "tri"){   #tmendelian error of third allele that cannot be consistent with any parent
             $genotype = 33; 
           } else { #biallelic error either sex/auto
            if ($sexLinkage eq "sex"){ #biallelic mendelian error sex chromosome
              if ($sex==$homogameticSexCode){$genotype = simHomogameticSexChrMendelianError($parent1,$parent2); 
              } else {
              $genotype = simHeterogameticSexChrMendelianError($parent1,$parent2); 
              }
            } else { # biallelic mendelian error on autosome. 
            $genotype = simMendelianError($parent1,$parent2); # this error mechanism introduces a mendelian error that is consistent with the biallelic gametes of the parents
            }
          } #end biallelic error
        }#end simulate error loop
          if (rand() < $missingness){ #set the genotypes for this individual to missing
           @genoVector = ($genoMissing,$genoMissing,$genoMissing,$genoMissing,$genoMissing,$genoMissing,$genoMissing,$genoMissing,$genoMissing,$genoMissing);
          } else {
            if ($genotype == 11){ $genoVector[0] = 0; $genoVector[1] =sprintf("%.0f",$expOffset+random_exponential(1,$lambda1)); $genoVector[4] = sprintf("%.0f",$expOffset+random_exponential(1,$lambda2));}
            if ($genotype == 12){ $genoVector[1] = 0; $genoVector[0] =sprintf("%.0f",$expOffset+random_exponential(1,$lambda1)); $genoVector[4] = sprintf("%.0f",$expOffset+random_exponential(1,$lambda1));}
            if ($genotype == 22){ $genoVector[4] = 0; $genoVector[1] =sprintf("%.0f",$expOffset+random_exponential(1,$lambda1)); $genoVector[0] = sprintf("%.0f",$expOffset+random_exponential(1,$lambda2));}
            if ($genotype == 33){ $genoVector[7] = 0; $genoVector[8] =sprintf("%.0f",$expOffset+random_exponential(1,$lambda1)); $genoVector[9] = sprintf("%.0f",$expOffset+random_exponential(1,$lambda2));}
          }
                                    #AA  AC  AG  AT  CC  CG  CT  GG  GT  TT
                                    # 0   1   2   3   4   5   6   7   8   9
} #end parent/offpring if/else

my $genoVector = join("\t",@genoVector);
my @ary = ();
if ($individual == 1){$ary[0] = "fam".$fam."FA"; $ary[1]=0;$ary[2]=0;$ary[3]=$homogameticSexCode;}
if ($individual == 2){$ary[0] = "fam".$fam."MO"; $ary[1]=0;$ary[2]=0;$ary[3]=$heterogameticSexCode;}
if (rand() < $sexUnknown){$sex = 0;}
if ($individual > 2){$ary[0] = "fam".$fam."OFF".$offspringNum; $ary[1]="fam".$fam."FA";$ary[2]="fam".$fam."MO";$ary[3]=$sex;}

#this output format is what we're after
#Chr	Pos	SNP	Fam	IndivID	Dad	Mom	Sex	AA	AC	AG	AT	CC	CG	CT	GG	GT	TT

  my $printHumanReadable = 0;

  if ($printHumanReadable == 1){
           # Output that Cris likes 
        if ($ary[1] =~ "fam"){
           print OUT "$chr\t1\t$snp\t$fam\t$ary[0]\t$ary[1]\t$ary[2]\t$ary[3]\t\t\t$genoVector\n";
        } else {
           print OUT "$chr\t1\t$snp\t$fam\t$ary[0]\t$ary[1]\t$ary[2]\t\t\t$ary[3]\t\t\t$genoVector\n";
        }
  } else {
         # Standard output format
         print OUT "$chr\t1\t$snp\t$fam\t$ary[0]\t$ary[1]\t$ary[2]\t$ary[3]\t$genoVector\n";
  }

} #end foreach individual in nuclear family


                 } #end foreach family, exiting bat territory
               } #end foreach iteration
             }#end foreach etype
            }#end foreach sex
          } #end foreach numFams
        } #end foreach missingness
      } #end foreach setCV
    } #end foreach numMendErrPerFam
  } #end foreach numOff
} #end foreach maf

} #end foreach option  #bazinga


################################################################################# subroutines and other stuff

sub simParent { #require a maf, return a founder genotype from the specified founder allele frequency
 if (scalar(@_ != 1)){die "please feed subroutine simParent a single minor allele frequency\n";}
 my $maf = $_[0];
 if ($maf > 0.5){die "please feed subroutine simParent a valid minor allele frequency\n";}
 my $allele1 = 1;
 my $allele2 = 1;
 if (rand()<$maf){$allele1 = 2;}
 if (rand()<$maf){$allele2 = 2;}
 my $geno = $allele1.$allele2;
 if ($allele1 > $allele2){$geno = $allele2.$allele1;}
 return($geno);
}

sub simHeterogameticParent { #require a maf, return a founder genotype from the specified founder allele frequency 
                    #this is not drawn from hardy weinberg expectation, but is the pseudo dipoid genotype of a truly heterogametic sex chromosome snp
 if (scalar(@_ != 1)){die "please feed subroutine simHomoParent a single minor allele frequency\n";}
 my $maf = $_[0];
 if ($maf > 0.5){die "please feed subroutine simHomoParent a valid minor allele frequency\n";}
 my $allele1 = 1;
 if (rand()<$maf){$allele1 = 2;}
 my $geno = $allele1.$allele1;
 return($geno);
}

sub simOffspring { #require genotypes for 2 parents, return a possible offspring genotype consistent mendelian transmission
 if (scalar(@_ != 2)){die "please feed subroutine simOffspring genotypes for two parents\n";}
 my @parent1 = split(//,$_[0]); if (scalar(@parent1 != 2)){die "please feed subroutine simOffspring biallelic genotypes for parent1\n";}
 my @parent2 = split(//,$_[1]); if (scalar(@parent2 != 2)){die "please feed subroutine simOffspring biallelic genotypes for parent2\n";}
 my $gamete1 =(); my $gamete2 =();
 if (rand()>.5){$gamete1 = $parent1[0]}else{$gamete1 = $parent1[1];}
 if (rand()>.5){$gamete2 = $parent2[0]}else{$gamete2 = $parent2[1];}
 my $offspring = $gamete1.$gamete2;
 if ($gamete1 > $gamete2){$offspring = $gamete2.$gamete1;};
 return($offspring);
}

sub simHeterogameticOffspring { #require genotypes for 2 parents, return a possible offspring genotype consistent mendelian transmission on a sex linked chromosome
                      #for the time being, assume parent 2 is heterogametic, and the offspring is heterogametic
 if (scalar(@_ != 2)){die "please feed subroutine simSexOffspring genotypes for two parents\n";}
 my @parent1 = split(//,$_[0]); if (scalar(@parent1 != 2)){die "please feed subroutine simSexOffspring biallelic genotypes for parent1\n";}
 my @parent2 = split(//,$_[1]); if (scalar(@parent2 != 2)){die "please feed subroutine simSexOffspring biallelic genotypes for parent2\n";}
 my $gamete1 =(); my $gamete2 =();
 if (rand()>.5){$gamete1 = $parent1[0]}else{$gamete1 = $parent1[1];}
 $gamete2 = $gamete1;  
 my $offspring = $gamete1.$gamete2;  #in truth, the offspring has only one allele of this diplo-genotype
 if ($gamete1 > $gamete2){$offspring = $gamete2.$gamete1;};
 return($offspring);
}

sub simMendelianError { #require genotypes for 2 parents, return an incompatible offspring genotype, otherwise die
 if (scalar(@_ != 2)){die "please feed subroutine simMendelianError genotypes for two parents\n";}
# my @parent1 = split(//,$_[0]); if (scalar(@parent1 != 2)){die "please feed subroutine simOffspring biallelic genotypes for parent1\n";}
# my @parent2 = split(//,$_[1]); if (scalar(@parent2 != 2)){die "please feed subroutine simOffspring biallelic genotypes for parent2\n";}
   my $parent1 = $_[0];
   my $parent2 = $_[1];
   my $errorGenotype = ();
   if ($parent1 eq "12" && $parent2 eq "12"){die "Dude. Check your math in the simMendelianError subroutine.  Two parents are hets.\n";}
   if ($parent1 eq "11" && $parent2 eq "11"){if (rand()<.5){$errorGenotype = "12";}else{$errorGenotype = "22";}}
   if ($parent1 eq "11" && $parent2 eq "12"){$errorGenotype = "22";}
   if ($parent1 eq "11" && $parent2 eq "22"){if (rand()<.5){$errorGenotype = "11";}else{$errorGenotype = "22";}}
   if ($parent1 eq "12" && $parent2 eq "11"){$errorGenotype = "22";}
   if ($parent1 eq "12" && $parent2 eq "22"){$errorGenotype = "11";}
   if ($parent1 eq "22" && $parent2 eq "11"){if (rand()<.5){$errorGenotype = "11";}else{$errorGenotype = "22";}}
   if ($parent1 eq "22" && $parent2 eq "12"){$errorGenotype = "11";}
   if ($parent1 eq "22" && $parent2 eq "22"){if (rand()<.5){$errorGenotype = "11";}else{$errorGenotype = "12";}}
 return($errorGenotype);
}

sub simHomogameticSexChrMendelianError { #require genotypes for 2 parents, return an incompatible offspring genotype, otherwise die
 if (scalar(@_ != 2)){die "please feed subroutine simHomogameticSexChrMendelianError genotypes for two parents\n";}
   my $parent1 = $_[0];
   my $parent2 = $_[1];
   my $errorGenotype = ();
   if ($parent2 eq "12"){die "Dude. Check your math in the simHomogameticSexChrMendelianError subroutine.  heterogametic parent is a het.\n";}
   if ($parent1 eq "11" && $parent2 eq "11"){if (rand()<.5){$errorGenotype = "12";}else{$errorGenotype = "22";}}
   if ($parent1 eq "12" && $parent2 eq "11"){$errorGenotype = "22";}
   if ($parent1 eq "22" && $parent2 eq "11"){if (rand()<.5){$errorGenotype = "11";}else{$errorGenotype = "22";}}
   if ($parent1 eq "11" && $parent2 eq "22"){if (rand()<.5){$errorGenotype = "11";}else{$errorGenotype = "22";}}
   if ($parent1 eq "12" && $parent2 eq "22"){$errorGenotype = "11";}
   if ($parent1 eq "22" && $parent2 eq "22"){if (rand()<.5){$errorGenotype = "11";}else{$errorGenotype = "12";}}
 return($errorGenotype);
}

sub simHeterogameticSexChrMendelianError { #require genotypes for 2 parents, return an incompatible offspring genotype, otherwise die
 if (scalar(@_ != 2)){die "please feed subroutine simHeterogameticSexChrMendelianError genotypes for two parents\n";}
   my $parent1 = $_[0];
   my $parent2 = $_[1];
   my $errorGenotype = ();
   if ($parent1 eq "12" or $parent2 eq "12"){die "Dude. Check your math in the simHeterogameticSexChrMendelianError subroutine.  a parent is a het.\n";}
   if ($parent1 eq "11"){$errorGenotype = "22";}
   if ($parent1 eq "22" ){$errorGenotype = "11";}
 return($errorGenotype);
}

#info on exponential distribution
#http://search.cpan.org/~grommel/Math-Random-0.70/Random.pm
#random_exponential($n, $av), When called in an array context, returns an array of $n deviates generated from the exponential distribution with mean $av. When called in a scalar context, generates and returns only one such deviate as a scalar, regardless of the value of $n.

sub gaussian_rand {  #from Perl cookbook,  polar Box Muller method = slow but easy
    my ($u1, $u2);  # uniformly distributed random numbers
    my $w;          # variance, then a weight
    my ($g1, $g2);  # gaussian-distributed numbers

    do {
        $u1 = 2 * rand() - 1;
        $u2 = 2 * rand() - 1;
        $w = $u1*$u1 + $u2*$u2;
    } while ( $w >= 1 );

    $w = sqrt( (-2 * log($w))  / $w );
    $g2 = sprintf("%.0f", ($u1 * $w));
    $g1 = sprintf("%.0f", ($u2 * $w));
    # return both if wanted, else just one
    return wantarray ? ($g1, $g2) : $g1;
}

sub truncZero { #truncate negative numbers to zero
  my $number = shift;
  if (0 > $number){ $number = 0;}
  return $number;
}

sub trim($){
	my $string = shift;
	$string =~ s/^\s+//;
	$string =~ s/\s+$//;
	return $string;
}

sub numerically {
    $a <=> $b;
}
sub round {my($number) = shift;
    return (int($number + .5 * ($number <=> 0)));
}

sub mean{
    my @ary = @_;
    my $n = scalar(@ary);
    my $sum = 0;
    foreach (@ary) {
	  if (!defined $_ || $_ eq ""){
	     $n--;
		 next;
	  }else{
	  $sum += $_ ;
	  }
    }
    if ($n == 0){ return ("");}
    if ($n == 1){ return ("@ary");}
    return ($sum/$n);
}
sub median{
    my @ary = ();
	my @ary2 = @_;
	foreach my $i (0 .. scalar(@ary2)-1){
	 if (defined $ary2[$i] && $ary2[$i] ne ""){
	    push(@ary,$ary2[$i]);
	 }
	}
    @ary = sort numerically @ary;
    my $n = scalar(@ary);
    if ($n == 0){ return ("");}
    if ($n == 1){ return ("@ary");}
    my $mid = int($n/2);
    my $median = ($n % 2) ? $ary[$mid] : ($ary[$mid] + $ary[$mid-1])/2;
    return $median;
}
sub max{
    my @ary = sort numerically @_;
    my $n = scalar(@ary);
    if ($n == 0){ return ("");}
    if ($n == 1){ return ("@ary");}
    return $ary[-1];
}
sub min{
    my @ary = sort numerically @_;
    my $n = scalar(@ary);
    if ($n == 0){ return ("");}
    if ($n == 1){ return ("@ary");}
    my $min = shift(@ary);
    return $min;
}
sub variance{
    my @ary2 = @_;
	my @ary = ();
	foreach my $a (@ary2){
	  if (defined $a && $a ne "" && $a ne " "){
  	    push(@ary,$a);
	  }
	}
    my $n = scalar(@ary);
    if ($n == 0){ return ("");}
    if ($n == 1){ return ("NA");}
    my $SSE = 0;
    my $mean=mean(@ary);
    foreach (@ary) {
        my $deviation = $_ - $mean;
        $SSE += $deviation**2;
#        if ($variance) {
#	    $skew /= ($n * $variance * $standard_deviation);
#	    $kurtosis = $kurtosis/($n * $variance * $variance) - 3.0;
#        }

    }
   return ($SSE/($n-1));
}
sub stdev{ return((variance(@_)^.5))}
#sub Q1{
#    my @ary = sort numerically @_;
#    my $n = scalar(@ary);
#    if ($n == 0){ return ("");}
#    if ($n => 4){ return ("NA");}
#    my $q25 = int($n/4);
#    my $Q25 = ($n % 2) ? $ary[$q25] : ($ary[$q25] + $ary[$q25-1])/2;
#    return $Q25;
#}

#sub Q3{
#    my @ary = sort numerically @_;
#    my $n = scalar(@ary);
#    if ($n == 0){ return ("");}
#    if ($n => 4){ return ("NA");}
#    my $q75 = int(3*$n/4);
#    my $Q75 = ($n % 2) ? $ary[$q75] : ($ary[$q75] + $ary[$q75-1])/2;
#    return $Q75;
#}
