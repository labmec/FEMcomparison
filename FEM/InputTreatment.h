//
// Created by victor on 02/09/2020.
//

#ifndef FEMCOMPARISON_INPUTTREATMENT_H
#define FEMCOMPARISON_INPUTTREATMENT_H

#include "DataStructure.h"


void EvaluateEntry(int argc, char *argv[],PreConfig &eData);
void ReadEntry(ProblemConfig &config, PreConfig &preConfig);
void InitializeOutstream(PreConfig &eData,char *argv[]);
void IsInteger(char *argv);
void Configure(ProblemConfig &config,int ndiv,PreConfig &pConfig,char *argv[]);
void CharReplace(std::string &str, char find, char replace );
void InitializeSpeedUp(PreConfig &pConfig);

#endif //FEMCOMPARISON_INPUTTREATMENT_H
