#include "../GridFunctions/grid_functions.h"
#include "../ElementData/d4est_element_data.h"
#include "../dGMath/d4est_operators.h"
#include "../LinearAlgebra/d4est_linalg.h"
#include "../Utilities/d4est_util.h"
#include "../Flux/d4est_mortars_compute_flux.h"
#include <curved_dg_norm.h>
#include <ip_flux_params.h>

/* See houston2007 for the definition of the dg-norm */
