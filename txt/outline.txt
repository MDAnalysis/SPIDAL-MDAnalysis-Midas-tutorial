.. -*- coding: utf-8 -*-

=========================================
 Tutorial for MDAnalysis + radical.pilot
=========================================

People
======

- Shantenu Jha
- Ioannis (Giannis) Paraskevakos i.paraskev@rutgers.edu
- Oliver
- Sean, Mahzad?

Geoffrey Fox will present the tutorial.


Existing work
=============

PSA
---
- implementation in MDA 0.16.0 dev or 0.15.0
- datasets
  - AdK dims/froda (big ensemble)
  - Diptheriatoxin dims/froda (big ensemble) (best example, I think)



RP
--

Giannis' scripts can be found here:
https://github.com/radical-cybertools/midas
 
Please, open new issues for everything that may be need either a quick
extension or more detailed explanation.

https://github.com/radical-cybertools/midas/tree/master/HausdorffDistance


Outline
=======

PARALLEL ANALYSIS OF MACROMOLECULAR TRANSITIONS WITH MDANALYSIS AND RADICAL.PILOT

The tutorial will showcase how to use radical.pilot to parallelize the analysis of an ensemble of macromolecular transition paths with Path Similarity Analysis (PSA).

BACKGROUND

- introduce problem: 
	- quantitative comparison of transitions: PSA
	- large ensembles and all-pairs: pleasingly parallel with radical.pilot (advantage: existing infrastructure, just run serial jobs)

- describe data
	- MD trajectories
	- numbers, sizes

- describe algorithms (?) – briefly
	- PSA: rmsd fitting and metric calculation (Hausdorff or Fréchet)
	- specifics of the partitioning of the ensemble for radical.pilot

- description of frame works
	- MDAnalysis (brief)
	- radical.pilot – or more likely you will have done this previously ?


TUTORIAL

- serial job
	- MDAnalysis / python script (explain basic functionality and relate to algorithm)
	- shows the basic work unit
	- makes it easier to understand how we then use rp to run many of them

- set up
	- data organization
	- radical.pilot script

- running rp

- analysis
	- gathering data
	- plot and interpretation
