#!/usr/bin/perl -w

use strict;
use warnings;

############################################################################
# mendelian error simulation for genotype vector data on nuclear pedigrees
# written by CVH, cvv8@cornell.edu beginning about Sep 19, 2012.
# v2: -simulate phred scores based on exponential distribution
#     -simulate multiple iterations for each parameter configuration
# v3: -change mendelian errors in offspring from 'non-biallelic' model to 'biallelic mendelian incompatible' model
# v4: -add capability to simulate sex linked snps March 1, 2013
# v5: -add mendelian error model to sex chromosomes March 7, 2013
# v6: -change the mendelian error model to 'biallelic' or 'triallelic' variable mendelian error types, also correct missingness to be independent of mendelian error status, i.e. with missingness, some pedigrees could be missing instead of having an error
# V7-8: -change the way that offspring genotypes are picked from parents, incorporate bivariate uniform sampler, allow for double het parents, disallow specification of the number of mendelian errors
# to do, a refactor should include subroutines for simulating genotypes(founder,het/homogametic,ParentalGenos) and simError(het/homogametic,sex,ParentalGenos), this moves the logic to the subroutines
############################################################################
#parameters and description

#these are user defined variables, it is ok to change these
my $sexRatio = .5;             #the probability the offspring will be heterogametic.
my $sexUnknown = 0;            #the probability that an offspring with masked (unkonwn) sex
my $heterogameticSexCode = 2;  #code for heterogametic sex, 2 for ZW systems, 1 for XY
my $homogameticSexCode = 1;    #code for homogametic sex, 1 for ZZ systems, 2 for XX, or whatever...
my $genoNull = "";             #null character for genotypes without phred scores, i.e. "" or "."(NOT missing data)
my $genoMissing = "-1";        #symbol for genotypes of individuals who are missing all genotypes
my $iterations = 1;         #run each of the parameter sets this number of times
my $printHumanReadable = 0;    #0 for standard output for MendelChecker, set to 1 for human friendly output (cannot be parsed by MendelChecker)



### global variables, do not change these
my @sexLinkage = ();           #sex linked or autosomal transmission
my @errorType = ();            #error model for mendelian error, either biallelic and inconsistent with mendelian transmission, or triallelic where the mendelian error is an offspring homozygous for the third allele
my @MAFS = ();                 #minor allele frequency of the putative biallelic SNP
my @numOffspring = ();         #number of offspring per family
my @numMendErrPerFam = ();     #number of mendelian errors per family, always in the children
my @setGeno = ();              #set the coefficient of variation (as a percent) with this parameter
my @missingness = ();          #zero to one, this proportion of individuals will have no data
my @numFams = ();              #this is the number of pedigrees to simulate for each condition
my @pOffIndep = ();            #probability that genotypes for 'offspring' are truly independent of parents
my $optionLabel = ();          #you guessed it, a label to collect all these variables
my $corr = 0;                  #correlation between parents and offspring, (keep at 0);
my $snp = 0;                   #simple counter
my $file = shift;              #output file header
   $file = "test";
#############################################################################

open (KEY, ">./$file.simFSJGenotypesv10.key") or die;
print KEY "SNP\tParameter\tParameterValue\tSex\tIteration\tParameterSet\n";
open (FAMKEY, ">./$file.simFSJGenotypesv10.famkey") or die;
print FAMKEY "SNP\tFAM\tParameter\tParameterValue\n";
open (OUT, ">./$file.simFSJGenotypesv10.txt") or die;
print OUT "#CHR\tPOS\tSNP\tFAM\tID\tFA\tMO\tSEX\tAA\tAC\tAG\tAT\tCC\tCG\tCT\tGG\tGT\tTT\n";

foreach my $option (9){ #iterate through all of the parameters of interest, or don't by changing this variable
#foreach my $option (1 .. 9){ #iterate through all of the parameters of interest, or don't by changing this variable
#parameter vectors, when the parameter in question is the random variable, this is range of values, else the default value
if ($option == 0){ $optionLabel = "NONE";}
if ($option == 1){ $optionLabel = "MAF"; @MAFS = ("0.01","0.05","0.1","0.25","0.5");} else {@MAFS = ("0.01" ,"0.02" ,"0.03" ,"0.04" ,"0.05" ,"0.06" ,"0.07" ,"0.08" ,"0.09" ,"0.1" ,"0.11" ,"0.12" ,"0.13" ,"0.14" ,"0.15" ,"0.16" ,"0.17" ,"0.18" ,"0.19" ,"0.2" ,"0.21" ,"0.22" ,"0.23" ,"0.24" ,"0.25" ,"0.26" ,"0.27" ,"0.28" ,"0.29" ,"0.3" ,"0.31" ,"0.32" ,"0.33" ,"0.34" ,"0.35" ,"0.36" ,"0.37" ,"0.38" ,"0.39" ,"0.4" ,"0.41" ,"0.42" ,"0.43" ,"0.44" ,"0.45" ,"0.46" ,"0.47" ,"0.48" ,"0.49" ,"0.5")}
#if ($option == 1){ $optionLabel = "MAF"; @MAFS = ("0.01","0.05","0.1","0.25","0.5");} else {@MAFS = ("0.01","0.05","0.1","0.25",".5");}
if ($option == 2){ $optionLabel = "OFF"; @numOffspring = ("1","2","5","10");} else {@numOffspring = ("1");} #be merciful, keep these parameters low for high iterations
if ($option == 3){ $optionLabel = "ERR"; @numMendErrPerFam = ("0","1","2");} else { @numMendErrPerFam = ("0");}
if ($option == 4){ $optionLabel = "GENO";  @setGeno = ("3000","500","100","1");} else {@setGeno = ("3000");}
if ($option == 5){ $optionLabel = "MIS"; @missingness = ("0","0.1","0.2","0.3","0.4","0.5","0.6","0.7","0.8","0.9","1");} else {@missingness = ("0");}
if ($option == 6){ $optionLabel = "FAM"; @numFams = ("1","2","3","4","5","10","25","50");} else {@numFams = ("25");}
if ($option == 7){ $optionLabel = "SEX"; @sexLinkage = ("auto","sex");} else {@sexLinkage = ("auto","sex");}
if ($option == 8){ $optionLabel = "ETYPE"; @errorType = ("bi","tri");} else {@errorType = ("bi");}
if ($option == 9){ $optionLabel = "PINDEP"; @pOffIndep = ("0","1");} else {@pOffIndep = ("0");}


foreach my $maf (@MAFS){
  foreach my $numOff (@numOffspring){
    foreach my $numMendErrPerFam (@numMendErrPerFam){
      foreach my $setGeno (@setGeno){
        foreach my $missingness (@missingness){
          foreach my $numFam (@numFams){
           foreach my $sexLinkage (@sexLinkage){
            foreach my $eType (@errorType){
             foreach my $pindep (@pOffIndep){

            if ($numMendErrPerFam > 0 && $pindep > 0){die "this set of parameters has not been tested.\n";}

             foreach my $iteration (1 .. $iterations){
             $snp++;

            my $chr = "chr_".$optionLabel."_iter".$iteration."_maf".$maf."_off".$numOff."_fam".$numFam."_err".$numMendErrPerFam."_geno".$setGeno."_mis".$missingness."_sex".$sexLinkage."_eType".$eType."_pindep".$pindep;
            print "option:$optionLabel iter:$iteration ParameterSet:$chr\n";
            my $paraValue = '';
            if ($optionLabel eq "NONE"){
            	print KEY "$snp\t$optionLabel\tnone\t$sexLinkage\t$iteration\t$chr\n";
            	$paraValue = "none";
            }
            if ($optionLabel eq "MAF"){
            	print KEY "$snp\t$optionLabel\t$maf\t$sexLinkage\t$iteration\t$chr\n";
            	$paraValue = $maf;
            }
            elsif ($optionLabel eq "OFF"){
            	print KEY "$snp\t$optionLabel\t$numOff\t$sexLinkage\t$iteration\t$chr\n";
            	$paraValue = $numOff;
            }
            elsif ($optionLabel eq "ERR"){
            	print KEY "$snp\t$optionLabel\t$numMendErrPerFam\t$sexLinkage\t$iteration\t$chr\n";
            	$paraValue = $numMendErrPerFam;
            }
            elsif ($optionLabel eq "GENO"){
            	print KEY "$snp\t$optionLabel\t$setGeno\t$sexLinkage\t$iteration\t$chr\n";
            	$paraValue = $setGeno;
            }
            elsif ($optionLabel eq "MIS"){
            	print KEY "$snp\t$optionLabel\t$missingness\t$sexLinkage\t$iteration\t$chr\n";
            	$paraValue = $missingness;
            }
            elsif ($optionLabel eq "SEX"){
            	print KEY "$snp\t$optionLabel\t$sexLinkage\t$sexLinkage\t$iteration\t$chr\n";
            	$paraValue = $sexLinkage;
            }
            elsif ($optionLabel eq "ETYPE"){
            	print KEY "$snp\t$optionLabel\t$sexLinkage\t$sexLinkage\t$iteration\t$chr\n";
            	$paraValue = $eType;
            }
            elsif ($optionLabel eq "PINDEP"){
            	print KEY "$snp\t$optionLabel\t$sexLinkage\t$sexLinkage\t$iteration\t$chr\n";
            	$paraValue = $pindep;
            }else{
            	print KEY "$snp\t$optionLabel\t$numFam\t$sexLinkage\t$iteration\t$chr\n";
            	$paraValue = $numFam;
            }

#set up phred score simulation parameterized by exponential random variable with mean values parameterized below
my $lambda1 = $setGeno/10;  #this factor 'default 10' sets the parameter for the exponential distribution for the penultimate most likely genotype as ten times lower (phred score is 10 times lower) than the parameter for the third most likely
my $lambda2 = ($setGeno);
my $expOffset = 0; #offset for the phred score distribution, most likely is 0, penultimate is 'expOffset + exp(lambda1)', third is is 'expOffset + exp(lambda2)'
if ($lambda1 == 300){$expOffset = 100;}  #this is cludgy, second offset for high quality genotype configuration
if ($lambda1 == 50){$expOffset = 50;}  #this is cludgy, second offset for med quality genotype configuration

foreach my $fam (1 .. $numFam){  #start the inner most simulation loop for the nuclear family..  bat territory
# print "label $optionLabel chr $chr fam $fam\n";

print FAMKEY "$snp\t$fam\t$optionLabel\t$paraValue\n";

my $famSize = $numOff+2;
my $parent1 = ();
my $shadow1 = ();
my $parent2 = ();
my $shadow2 = ();
my $sex = 0;
foreach my $individual (1 .. $famSize){
       my $genotype = ();
       my @genoVector = ($genoNull,$genoNull,$genoNull,$genoNull,$genoNull,$genoNull,$genoNull,$genoNull,$genoNull,$genoNull);
       my $offspringNum = $individual-2; #which offspring, First?  Second? this variable is ignored unless positive
  if ($individual <= 2){ #individual is a parent/founder
       if ($individual == 1){
          if ($sexLinkage eq "sex" && $numMendErrPerFam > 0){  #simulating sex linked errors, we must force parent1 to be homozygous
             my $passParentalHet = 0;
             while ($passParentalHet == 0){  #force parent two to be a non-het, but still simulate from the specified population allele frequency
                   my @temp = simParent($maf,$corr);
                   $parent1 = $temp[0];
                   $shadow1 = $temp[1];
                  if ($parent1 ne "12"){
                     $passParentalHet = 1;
                  }
             }
          } else {
          my @temp = simParent($maf,$corr);
          $parent1 = $temp[0];
          $shadow1 = $temp[1];
          }
        $genotype = $parent1;
       }

       if ($individual == 2){
          if ($sexLinkage eq "sex"){  #simulating sex linked snps, force parent2 to be homozygous. Not Hardy Weinberg proportions  In truth, this parent has only one allele of the diplo-genotype.
#          if ($numMendErrPerFam > 0){ die "mendelian errors on sex chromosomes has not been implemented\n";}
            my @temp = simHeterogameticParent($maf,$corr);
            $parent2 = $temp[0];
            $shadow2 = $temp[1];
          } else {
          my $passParentalHets = 1;
          if ($numMendErrPerFam > 0){ $passParentalHets = 0;}
              my @temp = simParent($maf,$corr);
              $parent2 = $temp[0];
              $shadow2 = $temp[1];
          if ($parent1 eq "12"){  #if parent one is a het...
             while ($passParentalHets == 0){  #if forcing mendelian errors, we must force parent two to be a non-het, but still simulate from the specified population allele frequency
                my @temp = simParent($maf,$corr);
                $parent2 = $temp[0];
                $shadow2 = $temp[1];
                  if ($parent2 ne "12"){
                     $passParentalHets = 1;
                  }
             }
          }
         }
          $genotype = $parent2;
       } #end if individual 2
       if (rand() < $missingness){ #set the genotypes for this parent to missing
         @genoVector = ($genoMissing,$genoMissing,$genoMissing,$genoMissing,$genoMissing,$genoMissing,$genoMissing,$genoMissing,$genoMissing,$genoMissing);
       } else {
       if ($genotype == 11){ $genoVector[0] = 0; $genoVector[1] =sprintf("%.0f",$expOffset+random_exponential(1,$lambda1)); $genoVector[4] = sprintf("%.0f",$expOffset+random_exponential(1,$lambda2));}
       if ($genotype == 12){ $genoVector[1] = 0; $genoVector[0] =sprintf("%.0f",$expOffset+random_exponential(1,$lambda1)); $genoVector[4] = sprintf("%.0f",$expOffset+random_exponential(1,$lambda1));}
       if ($genotype == 22){ $genoVector[4] = 0; $genoVector[1] =sprintf("%.0f",$expOffset+random_exponential(1,$lambda1)); $genoVector[0] = sprintf("%.0f",$expOffset+random_exponential(1,$lambda2));}
       }
  }else{ #individual is an OFFSPRING, possibly heterogametic or homogametic, possibly a mendelian error or not, possibly not an error but missing
         #at the moment, sex chromosomes can be missing, but not errors
        if (rand() < $sexRatio){$sex=$heterogameticSexCode;}else{$sex=$homogameticSexCode;} #what about UNKNOWN?
        if ($sexLinkage eq "sex"){
           #this is where to check against a parameterized sex ratio, and determine the sex of the offspring and chromosomal configuration
          if ($sex == $heterogameticSexCode){
          if (rand() < $pindep){
            $genotype = simHeterogameticOffspring($shadow1,$shadow2);
          }else{
          $genotype = simHeterogameticOffspring($parent1,$parent2);
          }
          }else{
          if (rand() < $pindep){
            $genotype = simOffspring($shadow1,$shadow2);
          }else{
            $genotype = simOffspring($parent1,$parent2);
          }
          }
        } else { #offspring, autosomal snp
          if (rand() < $pindep){
          $genotype = simOffspring($shadow1,$shadow2);
          }else{
          $genotype = simOffspring($parent1,$parent2);
          }
        }
        if ($offspringNum <= $numMendErrPerFam){ #simulate a mendelian error in this offspring
           @genoVector = ($genoNull,$genoNull,$genoNull,$genoNull,$genoNull,$genoNull,$genoNull,$genoNull,$genoNull,$genoNull);
           if ($eType eq "tri"){   #tmendelian error of third allele that cannot be consistent with any parent
             $genotype = 33; 
           } else { #biallelic error either sex/auto
            if ($sexLinkage eq "sex"){ #biallelic mendelian error sex chromosome
             if ($sex==$homogameticSexCode){
              if (rand() < $pindep){
                $genotype = simHomogameticSexChrMendelianError($shadow1,$shadow2);
              }else{
                $genotype = simHomogameticSexChrMendelianError($parent1,$parent2);
              }
             } else { #offspring is homogametic
              if (rand() < $pindep){
               $genotype = simHeterogameticSexChrMendelianError($shadow1,$shadow2); 
              }else{
               $genotype = simHeterogameticSexChrMendelianError($parent1,$parent2);
              }
             }
            } else { # biallelic mendelian error on autosome.
             if (rand() < $pindep){
              $genotype = simMendelianError($shadow1,$shadow2); # this error mechanism introduces a mendelian error that is consistent with the biallelic gametes of the parents
             }else{
              $genotype = simMendelianError($parent1,$parent2);
             }
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



  if ($printHumanReadable == 1){
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
              } #end foreach corr (correlation between offspring and parent genotypes)
             }#end foreach etype
            }#end foreach sex
          } #end foreach numFams
        } #end foreach missingness
      } #end foreach setCV
    } #end foreach numMendErrPerFam
  } #end foreach numOff
} #end foreach maf

} #end foreach option  #bazinga


################################################################################# subroutines 

sub simParent { #require a maf, return a founder genotype from the specified founder allele frequency
 if (scalar(@_ != 2)){die "please feed subroutine simParent a minor allele frequency and correlation coefficient\n";}
 my $maf = $_[0]; my $corr = $_[1];
 if ($maf > 0.5|| $maf < 0){die "please feed subroutine simParent a valid minor allele frequency\n";}
 my $allele1 = 1; my $shadow1 = 1;
 my $allele2 = 1; my $shadow2 = 1;
 my $geno = (); my $shadow = ();
 if ($corr ==1){
    my $UDev = rand();
    my @BUDev = ($UDev,$UDev);
    if ($BUDev[0]<$maf){$allele1 = 2;}; if ($BUDev[1]<$maf){$shadow1 = 2;}
    $UDev = rand();
    @BUDev = ($UDev,$UDev); 
    if ($BUDev[0]<$maf){$allele2 = 2;}; if ($BUDev[1]<$maf){$shadow2 = 2;}
    $geno = $allele1.$allele2;
    $shadow = $shadow1.$shadow2;
    if ($allele1 > $allele2){$geno = $allele2.$allele1;}
    if ($shadow1 > $shadow2){$shadow = $shadow2.$shadow1;}
 } else {
    my @BUDev = (rand(),rand());
    if ($BUDev[0]<$maf){$allele1 = 2;}; if ($BUDev[1]<$maf){$shadow1 = 2;}
    @BUDev = (rand(),rand());
    if ($BUDev[0]<$maf){$allele2 = 2;}; if ($BUDev[1]<$maf){$shadow2 = 2;}
    $geno = $allele1.$allele2;
    $shadow = $shadow1.$shadow2;
    if ($allele1 > $allele2){$geno = $allele2.$allele1;}
    if ($shadow1 > $shadow2){$shadow = $shadow2.$shadow1;}
 }
 return($geno,$shadow);
}

sub simHeterogameticParent { #require a maf, return a founder genotype from the specified founder allele frequency 
                    #this is not drawn from hardy weinberg expectation, but is the pseudo dipoid genotype of a truly heterogametic sex chromosome snp
 if (scalar(@_ != 2)){die "please feed subroutine simHeterogameticParent a single minor allele frequency and correlation coefficient\n";}
 my $maf = $_[0]; my $corr = $_[1];
 if ($maf > 0.5|| $maf < 0){die "please feed subroutine simParent a valid minor allele frequency\n";}
 if ($corr < 0 || $corr > 1){die "please feed subroutine simParent a valid correlation coefficient\n";}
 my $allele1 = 1; my $shadow1 = 1;
 my $geno = (); my $shadow = ();
  if ($corr == 1){
    my $UDev = rand();
    my @BUDev = ($UDev,$UDev);
    if ($BUDev[0]<$maf){$allele1 = 2;}; if ($BUDev[1]<$maf){$shadow1 = 2;}
    $geno = $allele1.$allele1;
    $shadow = $shadow1.$shadow1;  
  } else {
    my @BUDev = (rand(),rand());
    if ($BUDev[0]<$maf){$allele1 = 2;}; 
    $geno = $allele1.$allele1;
    $shadow = $shadow1.$shadow1;
 }
 return($geno,$shadow);
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
 my $offspring = $gamete1.$gamete2;  #in truth, the offspring has only one allele, here it is coded as a diploid genotype
 if ($gamete1 > $gamete2){$offspring = $gamete2.$gamete1;};
 return($offspring);
}

sub simMendelianError { #require genotypes for 2 parents, return an incompatible offspring genotype, otherwise die
 if (scalar(@_ != 2)){die "please feed subroutine simMendelianError genotypes for two parents\n";}
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

#http://search.cpan.org/~grommel/Math-Random-0.70/Random.pm
sub random_exponential { # Arguments: ($n,$av), defaults (1,1)
    return wantarray() ? (genexp(1)) : genexp(1)
	if scalar(@_) == 0; # default behavior if no arguments
    my($n, $av) = @_;
    $av = 1 unless defined($av); # default $av is 1
#    croak "$av = \$av < 0\nin random_exponential(\$n,\$av)" if ($av < 0);
    return genexp($av) unless wantarray();
    my $val;
    my @ans = (0) x $n;
    foreach $val (@ans) { $val = genexp($av); }
    return @ans;
}

sub genexp {
  return(-$_[0]*log(rand())); #log means 'natural log' in Perl
}

sub truncZero { #truncate negative numbers to zero
  my $number = shift;
  if (0 > $number){ $number = 0;}
  return $number;
}

sub trim{
	my $string = shift;
	$string =~ s/^\s+//;
	$string =~ s/\s+$//;
	return $string;
}

sub numerically {
    $a <=> $b;
}

sub round {
    my $number = shift;
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

sub stdev{
   return((variance(@_)^.5))
}

