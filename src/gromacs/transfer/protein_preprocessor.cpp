#include "gromacs/transfer/transfer.h"

t_atoms* protein_preprocessor(t_atoms *atoms, t_state *state_global)
{
    int      i, j, rescounter;
    t_atoms *atoms2;
    /* build copy of ct_atoms */
    //snew(atoms2,1);
    //init_t_atoms(atoms2, atoms->nr, FALSE);
    atoms2 = copy_t_atoms(atoms);

    /* take a different part in the memory for the copied atomnames */
    snew(atoms2->atomname, atoms->nr);
    for (i = 0; i < atoms->nr; i++)
    {
        snew(atoms2->atomname[i], 1);
    }
    for (i = 0; i < atoms->nr; i++)
    {
        snew(atoms2->atomname[i][0], 6);
    }
//  for (i=0; i<atoms->nr; i++)
//  atoms2->atomname[i] = atoms2->atomname[0] + 6*i;
//  for (i=0; i<atoms->nr; i++)
//  sfree(atoms2->atomname[i]);
//  for (i=0; i<atoms->nr; i++)
//  snew(atoms2->atomname[i], atoms->nr);
//  for (i=0; i<atoms->nr; i++)
//  snew(atoms2->atomname[i][0], 6);
//atoms2->atomname[0][0]=&atoms2->atomname[0][0][0];

    printf("pointer %p   %p\n", atoms->atomname[2], atoms2->atomname[2]);
    printf("pointer %p   %p\n", atoms->atomname, atoms2->atomname);
    printf("pointer %p   %p\n", atoms->atomname[0], atoms2->atomname[0]);
    printf("pointer %p   %p\n", *(atoms->atomname[1]), *(atoms2->atomname[1]));
    printf("pointer %p   %p\n", atoms->atomname[1][0], atoms2->atomname[1][0]);
    printf("pointer %p\n", atoms2->atomname[0]);
    printf("pointer %p\n", *(atoms2->atomname[2]));

/*
    for (j=0; j<10; j++){
   printf("1name %s\n",*(atoms->atomname[j]));
   printf("1name %c %c %c %c\n",(atoms->atomname[j])[0][0], (atoms->atomname[j])[0][1],(atoms->atomname[j])[0][2],(atoms->atomname[j])[0][3]);
   printf("1name %c %c %c %c\n",(*(atoms->atomname[j]))[0], (*(atoms->atomname[j]))[1],(*(atoms->atomname[j]))[2],(*(atoms->atomname[j]))[3]);
   }
 */
    rescounter = -1;
    for (i = 0; i < atoms->nres; i++)
    {
        rescounter = rescounter+2;
        for (j = 0; j < atoms->nr; j++)
        {

//printf("resind %d\n",atoms->atom[j].resind);
//printf("1name %c %c %c %c\n",(atoms->atomname[j])[0][0], (atoms->atomname[j])[0][1],(atoms->atomname[j])[0][2],(atoms->atomname[j])[0][3]);
//printf("pointer %p   %p\n",*(atoms->atomname[j]),*(atoms2->atomname[j]));
//printf("1name %c %c %c %c\n",(*(atoms->atomname[j]))[0], (*(atoms->atomname[j]))[1],(*(atoms->atomname[j]))[2],(*(atoms->atomname[j]))[3]);
//printf("1name %s\n",*(atoms->atomname[j]));
            if (atoms->resinfo[atoms->atom[j].resind].nr-1 == i) /* atom j is in residue i. resind is index in resinfo starting from 0 .nr is number of residue starting from 1 */
            {
                if (!strcmp((*(atoms->atomname[j])), "C") ||
                    !strcmp((*(atoms->atomname[j])), "O")) // C and O
                {
                 // (*(atoms2->atomname[j])) = !strcmp((*(atoms->atomname[j])), "C") ? "CQMB" : "OQM";
                    if (!strcmp((*(atoms->atomname[j])), "C")) {
                      strcpy((*(atoms2->atomname[j])), "CQMB");
                    } else { 
                      strcpy((*(atoms2->atomname[j])), "OQM");
                    }
                    atoms2->atom[j].resind   = rescounter+1;
                }
                else if (!strcmp((*(atoms->atomname[j])), "CA") ||
                         !strcmp((*(atoms->atomname[j])), "HA"))
                {
                    atoms2->atom[j].resind   = rescounter-1;
                 // (*(atoms2->atomname[j])) = (!strcmp((*(atoms->atomname[j])), "C")) ? "CQMB" : "HQM";
                    if (!strcmp((*(atoms->atomname[j])), "C")) {
                      strcpy((*(atoms2->atomname[j])), "CQMB");
                    } else { 
                      strcpy((*(atoms2->atomname[j])), "HQM");
                    }
                    if (!strcmp((*(atoms->resinfo[atoms->atom[j].resind].name)), "PRO"))
                    {
                        atoms2->atom[j].resind = rescounter+1;
                    }                                                                                                  //merge res 0 1 2 for prolin
                }
                else if (!strcmp((*(atoms->atomname[j])), "N") ||
                         !strcmp((*(atoms->atomname[j])), "H"))
                {
                 // (*(atoms2->atomname[j])) = (!strcmp((*(atoms->atomname[j])), "N")) ? "NQMB" : "HQM";
                    if (!strcmp((*(atoms->atomname[j])), "N")) {
                      strcpy((*(atoms2->atomname[j])), "NQMB");
                    } else { 
                      strcpy((*(atoms2->atomname[j])), "HQM");
                    }
                    atoms2->atom[j].resind   = rescounter-1;
                }
                else
                {
                    printf("sidechain\n");
                    if (!strcmp((*(atoms->atomname[j])), "CG"))
                    {
                     // (*(atoms2->atomname[j])) = "CQMB";
                        strcpy((*(atoms2->atomname[j])), "CQMB");
                    }
                    else
                    {
//(*(atoms2->atomname[j]))=" ";
                        strncpy((*(atoms2->atomname[j])), (*(atoms->atomname[j])), 1);
                        strcat((*(atoms2->atomname[j])), "QM");
                        // (atoms2->atomname[j])[0][1]='Q'; //crop string. first letter is element type
                        //  (atoms2->atomname[j])[0][2]='M'; //crop string. first letter is element type
                        //  (atoms2->atomname[j])[0][3]='\0'; //crop string. first letter is element type
                        printf("name 2 xname %s  %s\n", (*(atoms->atomname[j])), (*(atoms2->atomname[j])));
                        //strcat((*(atoms2->atomname[j])), "QM");
                    }
                    atoms2->atom[j].resind = rescounter-1;
                    if (!strcmp((*(atoms->resinfo[atoms->atom[j].resind].name)), "PRO"))
                    {
                        atoms2->atom[j].resind = rescounter-3;
                    }                                                                                                  //merge res 0 1 2 for prolin
                }
            }
        }
    }
    atoms2->nres = rescounter;

    for (j = 0; j < atoms->nr; j++)
    {
        printf("final name %s %d\n", *(atoms2->atomname[j]), atoms2->atom[j].resind);
//printf("1name %c %c %c %c\n",(atoms->atomname[j])[0][0], (atoms->atomname[j])[0][1],(atoms->atomname[j])[0][2],(atoms->atomname[j])[0][3]);
//printf("1name %c %c %c %c\n",(*(atoms->atomname[j]))[0], (*(atoms->atomname[j]))[1],(*(atoms->atomname[j]))[2],(*(atoms->atomname[j]))[3]);
    }
    atoms = copy_t_atoms(atoms2);
    for (j = 0; j < atoms->nr; j++)
    {
        printf("final name %s %d\n", *(atoms->atomname[j]), atoms->atom[j].resind);
    }

    //write_sto_conf("preprocessor.gro", ".gro file to check residues",  atoms,  state_global->x[], NULL, NULL, NULL);
    FILE *f_ct_preprocessor = NULL;
    f_ct_preprocessor = fopen("preprocessor.gro", "w");
    fprintf(f_ct_preprocessor, "gro file with new residues\n %d\n", atoms->nr);
    for (j = 0; j < atoms->nr; j++)
    {
        fprintf(f_ct_preprocessor, "%5d%-5.5s%5.5s%5d %lf %lf %lf\n", atoms->atom[j].resind%100000, "NEW", *(atoms->atomname[j]), (j+1)%100000, state_global->x[j][XX], state_global->x[j][YY], state_global->x[j][ZZ]);
    }
    fprintf(f_ct_preprocessor, "10 10 10"); //pseudo-box
    //write_hconf_box(f_ct_preprocessor, -1, state_global->box);
    fclose(f_ct_preprocessor);

    return atoms;
}

