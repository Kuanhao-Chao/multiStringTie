#ifndef __PROCESSBUNDLE_C_H__
#define __PROCESSBUNDLE_C_H__

#pragma once

#include "global_params_C.h"
#include "infer_tranx_C.h"

void processBundle_CREATE_UNISPG(BundleData* bundle, UnispgGp_CREATE* unispg_gp, int fidx);

void noMoreBundles();

#endif