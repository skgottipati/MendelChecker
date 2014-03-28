###MendelChecker
a C++ software for quality control in next generation sequencing using Mendelian inheritance in pedigrees.
author: S. Gottipati <skgottipati@gmail.com>

---

###Requirements

GCC compiler version >4.7 is required.

### Compiling

Using wget to download the zip file on Unix/Linux machines

	$ wget https://github.com/skgottipati/MendelChecker/archive/master.zip
	$ unzip -a master.zip
	$ mv MendelChecker-master MendelChecker
	$ cd MendelChecker
	$ make

You could also checkout directly from github source repository for the latest version of the software.

	$ git clone git://github.com/skgottipati/MendelChecker.git
	$ cd MendelChecker
	$ make

###Command line options

Options:
	-h, --help            show this help message and exit
	--version             show program's version number and exit
	-f FILE, --genoped=FILE	input geno-ped file name with path
	-g FILE, --vcf=FILE	input vcf file name with path
	-e FILE, --ped=FILE	input ped file name with path
	-n INT, --snpsperloop=INT	number of snps compute per loop
	-d STRING, --genofield=STRING	VCF genotype field, options: PL, GL, GP (default: PL)
	-p DOUBLE, --sexPrior=DOUBLE	default: 0.05 sexPrior
	-u STRING, --uniform=STRING	default: false (population), true (uniform)

###Usage information


For more information, please visit http://code.google.com/p/mendelchecker/wiki/Documentation

###If you use MendelChecker, please cite 


    @misc{Van14MC,
      Author = {Van Hout, C.V., Chen N., Gotipati S., Clark A.G., },
      Title = { MendelChecker: Quality Control for Next Generation Sequencing Using Mendelian Inheritance in Pedigrees},
      Year  = {2014},
      Howpublished = {Bioinformatics}
    }

	@misc{Chen14MC,
      Author = {Chen N., Van Hout, C.V., Gotipati S., Fitzpatrick J., Clark A.G., },
      Title = { TBD},
      Year  = {2014},
      Howpublished = {TBD}
    }