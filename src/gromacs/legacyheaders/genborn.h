/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2008, The GROMACS development team.
 * Copyright (c) 2013, by the GROMACS development team, led by
 * Mark Abraham, David van der Spoel, Berk Hess, and Erik Lindahl,
 * and including many others, as listed in the AUTHORS file in the
 * top-level source directory and at http://www.gromacs.org.
 *
 * GROMACS is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public License
 * as published by the Free Software Foundation; either version 2.1
 * of the License, or (at your option) any later version.
 *
 * GROMACS is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with GROMACS; if not, see
 * http://www.gnu.org/licenses, or write to the Free Software Foundation,
 * Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA.
 *
 * If you want to redistribute modifications to GROMACS, please
 * consider that scientific software is very special. Version
 * control is crucial - bugs must be traceable. We will be happy to
 * consider code for inclusion in the official distribution, but
 * derived work must not be called official GROMACS. Details are found
 * in the README & COPYING files - if they are missing, get the
 * official version at http://www.gromacs.org.
 *
 * To help us fund GROMACS development, we humbly ask that you cite
 * the research papers on the package. Check out http://www.gromacs.org.
 */


#ifndef _genborn_h
#define _genborn_h

#include "typedefs.h"
#include "types/commrec.h"
#include "grompp.h"

#ifdef __cplusplus
extern "C" {
#endif

/* Still parameters - make sure to edit in genborn_sse.c too if you change these! */
#define STILL_P1  0.073*0.1              /* length        */
#define STILL_P2  0.921*0.1*CAL2JOULE    /* energy*length */
#define STILL_P3  6.211*0.1*CAL2JOULE    /* energy*length */
#define STILL_P4  15.236*0.1*CAL2JOULE
#define STILL_P5  1.254

#define STILL_P5INV (1.0/STILL_P5)
#define STILL_PIP5  (M_PI*STILL_P5)


/* Initialise GB stuff */
int init_gb(gmx_genborn_t **p_born,
            const t_commrec *cr, t_forcerec *fr, const t_inputrec *ir,
            const gmx_mtop_t *mtop, int gb_algorithm);


/* Born radii calculations, both with and without SSE acceleration */
int calc_gb_rad(t_commrec *cr, t_forcerec *fr, t_inputrec *ir, gmx_localtop_t *top, rvec x[], t_nblist *nl, gmx_genborn_t *born, t_mdatoms *md, t_nrnb     *nrnb);



/* Bonded GB interactions */
real gb_bonds_tab(rvec x[], rvec f[], rvec fshift[], real *charge, real *p_gbtabscale,
                  real *invsqrta, real *dvda, real *GBtab, t_idef *idef, real epsilon_r,
                  real gb_epsilon_solvent, real facel, const t_pbc *pbc, const t_graph *graph);



void gb_pd_send(t_commrec *cr, real *send_data, int nr);


/* Functions for setting up the F_GB12,13,14 lists in grompp */
int
init_gb_plist(t_params *p_list);

int
convert_gb_params(gmx_ffparams_t *ffparams, t_functype ftype, t_params *gb_plist, t_ilist *il);



/* Functions for calculating adjustments due to ie chain rule terms */
void
calc_gb_forces(t_commrec *cr, t_mdatoms *md, gmx_genborn_t *born, gmx_localtop_t *top,
               rvec x[], rvec f[], t_forcerec *fr, t_idef *idef, int gb_algorithm, int sa_algorithm, t_nrnb *nrnb,
               const t_pbc *pbc, const t_graph *graph, gmx_enerdata_t *enerd);


int
make_gb_nblist(t_commrec *cr, int gb_algorithm,
               rvec x[], matrix box,
               t_forcerec *fr, t_idef *idef, t_graph *graph, gmx_genborn_t *born);

void
make_local_gb(const t_commrec *cr, gmx_genborn_t *born, int gb_algorithm);

#ifdef __cplusplus
}
#endif

#endif /* _genborn_h */
