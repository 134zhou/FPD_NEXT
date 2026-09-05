#ifndef GLOBAL_DEFINITIONS_H
#define GLOBAL_DEFINITIONS_H

#define NUMBER_THREADS 4

/* type definitions */
typedef unsigned long long int steps;


const double DT = 0.5E-03;

const double RATIO_ETA = 50.0;

const double XI = 1.;
const double RADIUS = 3.2;

const double CELL_SIZE = 32.;
const double CUTOFF_POTENTIAL = 3.;

const int INTERVAL_POSITIONS = 1;
//const float INTERVAL_POSITIONS = 0.01;
const int INTERVAL_VELOSITY = 100;
//const int INTERVAL_VELOSITY = 1;

//const double g = -1.;//gravity

const double De = 50.0;
const double a = 1;
const double Re = 7.4;

const double sgm = 7.4;
const double eps = 57.1428571429;

const double rho = 7.4;
#ifdef _8BLOCKLOOP_
const int RANGE_CORE = 7;
#endif

#endif

