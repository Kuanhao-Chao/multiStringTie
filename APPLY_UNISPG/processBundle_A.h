#ifndef __PROCESSBUNDLE_A_H__
#define __PROCESSBUNDLE_A_H__

#pragma once

#include "global_params_A.h"
#include "infer_tranx_A.h"

void processBundle_APPLY_UNISPG(BundleData* bundle, UnispgGp_APPLY* unispgs);

void noMoreBundles();

#endif