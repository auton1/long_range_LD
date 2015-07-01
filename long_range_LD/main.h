/*
 * main.h
 *
 *  Created on: Feb 13, 2015
 *      Author: aauton
 */

#ifndef MAIN_H_
#define MAIN_H_

#include <cmath>
#include <fstream>
#include <iostream>
#include <limits>
#include <queue>
#include <random>
#include <sstream>
#include <string>
#include <tuple>
#include <vector>

#include "output_log.h"
#include "parameters.h"
#include "bcf_file.h"
#include "vcf_file.h"
#include "variant_file.h"
#include "header.h"

using namespace std;


long mod(long a, long b)
{ return (a%b+b)%b; }


#endif /* MAIN_H_ */
