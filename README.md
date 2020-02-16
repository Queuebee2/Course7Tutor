# Course7Tutor

### Plan (so far)
- create 'indices' from reference genome with `bowtie2-build`
- create alignment map with `bowtie2` from read date with the indices (creates SAM file)
- use samtools to convert the SAM to a BAM file.
- use samtools to sort, index and 'depth'. 

# Installing tools

### Java 11

#### Fix apt-get
`sudo apt-get update`

#### install java 11
`sudo apt install openjdk-11-jre-headless --fix-missing`
### IGV Browser

#### Download IGV
get it [here](https://software.broadinstitute.org/software/igv/download), download the _IGV and igvtools to run on the command line (all platforms)_

#### Install IGV
just unzip it
```unzip IGV_2.8.0.zip```
### Bowtie2 and Samtools
should be available from the advanced package tools thingy

```bash
sudo apt install bowtie2
sudo apt install samtools
```

#### Running IGV
Find the proper command in the readme file  
`java --module-path=lib -Xmx4g @igv.args --module=org.igv/org.broad.igv.ui.Main`  
is the way we did it.             

### manual install samtools
I tried manually installing samtools, it took me all this crap to get there.
```bash
sudo apt-get install libncurses5-dev libncursesw5-dev
sudo apt-get install libz-dev

sudo apt-get install python-dev
sudo apt-get install python-bzutils
sudo apt-get install libbz2-dev

sudo apt-get install -y liblzma-dev
sudo apt-get install libcurl4-openssl-dev

wget https://sourceforge.net/projects/samtools/files/samtools/1.10/samtools-1.10.tar.bz2/download
tar xvjf samtools-1.10.tar.bz2
cd samtools-1.10/
make
```
