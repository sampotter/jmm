#pragma once

#include "def.h"

bool face_in_cell(uint3 const f, uint4 const c);
bool point_in_face(size_t l, uint3 const f);
bool point_in_cell(size_t l, uint4 const c);
bool edge_in_face(uint2 const le, uint3 const lf);
int edge_cmp(uint2 const e1, uint2 const e2);
void R_from_n(dbl3 const n, dbl33 R);
