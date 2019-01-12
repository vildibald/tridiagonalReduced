#pragma once

// Switch between 32 bit (float) and 64 bit precision (double).
using number = float;
// using number = double;

// After switching to double, change functions sdttrfb and sdttrsb in utils.cpp to
// to ddttrfb and ddttrsb