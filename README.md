# Coalescence_Project
# TODO for this README: 
1. fix up the README 
2. read more of the papers and articles to better inform this README 
3. Fixing up neighbor Joining Tree

# Purpose
At the moment, this is a project to produce an coalescent tree using mosquito genomes to understand why mosquito are very adaptive when thy are actually terrible at adapting. I would very much like to understand how a seemingly structural variant led to such a huge dominio effect that led them becoming the biggest malaria vector.

# Background Information 
Based on Lobo and other papers, mosquito were mainly set to environments with a pretty wet climate or weather. Mosquitos were able to effectively dominate other types of enviornmnet due to their structural variants.

# What am I doing
There is a wide variety of methods I am employing to make this project work
1. Clustering algorithm based on similar information (gathered using BWA and Samtools)
2. Weighted phylogenetic tree or Neighbor Joining Tree (checks for the nearest distance and cluster from there)

This project folder will be based on neighbor joining tree.

# Steps to reproduce what I did
1. In progress at the moment as I'm still programming the code and what not
    1. More coming soon 

# Quick Q&A
1. Why not use a library for all your needs as python and R carry their best version of neighbor joining tree? 
    1. Dude, I'm just trying to learn more about evolutionary genomes 
    2. I do want to get better at programming so this will be a fun activity to learn how dynamic programming effectively makes this algorithm much better than O(n^5) (which is horrendous) 
2. Editors used or applications? 
    1. Just vscode (visual studio code)
        1. I'm pretty sure you could find all of my vscode setting in the hidden folder here to be honest
    2. UCSC hummingbird for computational support (they use slurm as their scheduling queue) 
    3. Sources of Papers are down at the citation below 
3. Is this solo or group? 
    1. Solo




# Citation or Sources that really helped (min. of 15 sources to be safe)
Disorganized but these are the ones I have been using to help me with this project 
1. The Lobo dude – the real dude 
2. Russ’s two papers – clustering paper and other inversino breakpoints
3. the other two mosquitos – background on how this works  
4. Coalescent papers (three probably to make sure I get it) 
5. GATK paper 
6. BLAST paper and how it works and why I should use what I use 
7. PICARD paper and how it works as well to understand 
8. BWA to understand why use mem and what not all the fun stuff,.
9. Neighborhood joining tree and bed files 
10. bed files
11. to 15. More information on neighbor joining tree and why they could be used at O(n^3) instead of O(n^5) 

# Note
1. There are three files in this directory: 
    1. numpy_nj.py
    2. distancetree.py
    3. distance_matrix.txt 
2. They are other files from other people just a sanity check for my sake as I attempt to build the algorithm from scratch and understand what might be wrong
    1. Also, I can't run the files if anyone is wondering and would require so much work to make it work which I might as well learn how to build the algorithm myself.