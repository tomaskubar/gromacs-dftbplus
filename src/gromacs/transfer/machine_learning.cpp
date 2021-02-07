/*************************
 * Neural Nets: Mila Kraemer
 *************************/

#include "gromacs/transfer/transfer.h"
//#include <Python.h>


// KERNEL RIDGE
void get_ml_hamiltonian(charge_transfer_t *ct, t_state *state_global, t_mdatoms *mdatoms);
void get_ml_hamiltonian(charge_transfer_t *ct, t_state *state_global, t_mdatoms *mdatoms)
{
    long    iatom, isite, jsite, count, i, j, k, natoms, ncm;
    double  *cm, *atomcharge, **pair_distance_matrix, **coordinates;
    double   gauss_kernel, norm2;

//  store the coordinates of atoms in QM zone
    for (iatom = 0; iatom < mdatoms->homenr; iatom++)
    {
        for (i = 0; i < ct->sites; i++)
        {
            for (j = 0; j < ct->site[i].atoms; j++)
            {
                if (iatom == ct->site[i].atom[j])
                {
                   for (k = 0; k < DIM; k++)
                       ct->coords[i][j][k] = state_global->x.rvec_array()[iatom][k]*10.0;
                }
            }
        }
    }

// calculate hamiltonian by fitting to a linear function by gaussian kernel regression
    natoms = 2*ct->site[0].atoms;
    ncm = natoms*(natoms+1)/2;

    snew(cm, ncm);
    snew(atomcharge, natoms);

    snew(pair_distance_matrix, natoms);
    for (i = 0; i < natoms; i++)
        snew(pair_distance_matrix[i], natoms);

    snew(coordinates, natoms);
    for (i = 0; i < natoms; i++)
        snew(coordinates[i], DIM);

// NOTE ASSUME ALL SITES ARE THE SAME!!
    for (i = 0; i < ct->site[0].atoms; i++)
    {
        switch (ct->site[0].atomtype[i])
        {
          case 0: atomcharge[i] = 6.0; break; // C
          case 1: atomcharge[i] = 1.0; break; // H
          case 2: atomcharge[i] = 7.0; break; // N
          case 3: atomcharge[i] = 8.0; break; // O
        }
        atomcharge[ct->site[0].atoms+i] = atomcharge[i];
    }

// compute the diagonal terms of pair distance matrix
    for (i = 0; i < natoms; i++)
        pair_distance_matrix[i][i] = 0.5 * pow(atomcharge[i], 2.4);

    for (isite = 0; isite < ct->sites; isite++)
        for (jsite = 0; jsite < ct->sites; jsite++)
            ct->hamiltonian_ml[isite][jsite] = 0.0;

// diagonal terms predictions
    natoms = ct->site[0].atoms;
    ncm = natoms*(natoms+1)/2;

    for (isite = 0; isite < ct->sites; isite++)
    {
        count = 0;

        for (i = 0; i < natoms; i++)
        {
            for (j = i + 1; j < natoms; j++)
            {
                norm2 = 0.0;
                for (k = 0; k < DIM; k++)
                    norm2 += SQR(ct->coords[isite][j][k] - ct->coords[isite][i][k]);

// off-diagonal terms pair distance matrix
                pair_distance_matrix[i][j] = atomcharge[i] * atomcharge[j] / sqrt(norm2);
                pair_distance_matrix[j][i] = pair_distance_matrix[i][j];
            }

// coulomb matrix
            for (j = 0; j <= i; j++)
                cm[count+j] = pair_distance_matrix[i][j];

            count += i+1;
        }

// gaussian kernel function
        for (i = 0; i < ct->ntrain[0]; i++)
        {
            norm2 = 0.0;

            for (j = 0; j < ncm; j++)
            {
                norm2 += SQR(cm[j] - ct->cm_train[0][i][j]);
            }

            gauss_kernel = exp(-norm2/2.0/SQR(ct->sigma[0]));

            ct->hamiltonian_ml[isite][isite] += ct->alpha[0][i]*gauss_kernel;
        }

// Make the predictions
        ct->hamiltonian_ml[isite][isite] /= HARTREE_TO_EV;

// implicit reorganization
        if (ct->jobtype == cteNOMOVEMENT)
        {
        ct->hamiltonian_ml[isite][isite] -= ct->lambda*ct->occupation[isite];
    }
        else
        {
            ct->hamiltonian_ml[isite][isite] -= ct->lambda*(SQR(ct->tfs_diab[isite]) + SQR(ct->tfs_diab[isite+ct->dim]));
        }

    }

// off-diagnoal terms
    natoms = 2*ct->site[0].atoms;
    ncm = natoms*(natoms+1)/2;

    for (isite = 0; isite < ct->sites - 1; isite++)
    {
        for (jsite = isite + 1; jsite < ct->sites; jsite++)
        {
            for (i = 0 ; i < ct->site[0].atoms; i++)
            {
                for (k = 0 ; k < DIM; k++)
                {
                    coordinates[i][k] = ct->coords[isite][i][k];
                    coordinates[i+ct->site[0].atoms][k] = ct->coords[jsite][i][k];

                }

            }

            count = 0;
            for (i = 0; i < natoms; i++)
            {
                for (j = i + 1; j < natoms; j++)
                {
                    norm2 = 0.0;

                    for (k = 0; k < DIM; k++)
                        norm2 += SQR(coordinates[j][k] - coordinates[i][k]);

// off-diagonal terms pair distance matrix
                    pair_distance_matrix[i][j] = atomcharge[i] * atomcharge[j] / sqrt(norm2);
                    pair_distance_matrix[j][i] = pair_distance_matrix[i][j];
                }

// coulomb matrix
                for (j = 0; j <= i; j++)
                    cm[count+j] = pair_distance_matrix[i][j];

                count += i+1;
            }

// gaussian kernel function
            for (i = 0; i < ct->ntrain[1]; i++)
            {
                norm2 = 0.0;

                for (j = 0; j < ncm; j++)
                    norm2 += SQR(cm[j] - ct->cm_train[1][i][j]);

                gauss_kernel = exp(-norm2/2.0/SQR(ct->sigma[1]));

                ct->hamiltonian_ml[isite][jsite] += ct->alpha[1][i]*gauss_kernel;
            }

// Make the predictions
            ct->hamiltonian_ml[isite][jsite] /= HARTREE_TO_EV;
            ct->hamiltonian_ml[jsite][isite] = ct->hamiltonian_ml[isite][jsite];
        }
    }

    for (i = 0; i < natoms; i++)
    {
        sfree(pair_distance_matrix[i]);
        sfree(coordinates[i]);
    }

    sfree(cm);
    sfree(atomcharge);
    sfree(pair_distance_matrix);
}


/* MACHINE LEARNING WITH TF NNs -- MK 2020 */

#if (GMX_TENSORFLOW)

/*TF tensors need this*/
void NoOpDeallocator(void* data, size_t a, void* b) {(void) data; (void) a; (void) b;}

void initialize_nn_model(tf_model *model, const char *path, const char *hyperpath) {
  /* Initialize saved model from folder path into struct. */
  model->graph = TF_NewGraph();
  model->status = TF_NewStatus();
  TF_SessionOptions* SessionOpts = TF_NewSessionOptions();
  TF_Buffer* RunOpts = NULL;
  const char* tags = "serve";
  model->session = TF_LoadSessionFromSavedModel(SessionOpts, RunOpts, path, &tags, 1, model->graph, NULL, model->status);
  printf("session loaded successfully!\n");
  TF_DeleteSessionOptions(SessionOpts);
  TF_DeleteBuffer(RunOpts);

  /* Preparation of Input and output tensors */
  model->num_inputs = 1; //one input tensor
  model->input = (TF_Output*) malloc(sizeof(TF_Output) * model->num_inputs);
  // serving_default_input_1 is the default name, get from saved_model_cli
  TF_Output t0 = {TF_GraphOperationByName(model->graph, "serving_default_input_1"), 0};
  model->input[0] = t0;

  // we have two output tensors - the coupling and the force
  model->num_outputs = 2;
  model->output = (TF_Output*) malloc(sizeof(TF_Output) * model->num_outputs);
  // StatefulPartitionedCall is the default name, get from saved_model_cli
  TF_Output t2 = {TF_GraphOperationByName(model->graph, "StatefulPartitionedCall"), 0};
  model->output[0] = t2;
  TF_Output t3 = {TF_GraphOperationByName(model->graph, "StatefulPartitionedCall"), 1};
  model->output[1] = t3;

  model->input_val = (TF_Tensor**)malloc(sizeof(TF_Tensor*)*model->num_inputs);
  model->output_val = (TF_Tensor**)malloc(sizeof(TF_Tensor*)*model->num_outputs);

  FILE *hyperfile = NULL;
  hyperfile = fopen(hyperpath, "r");
  float hypers[5];
  printf("Reading hyperparameter file\n");
  for (int i=0; i<5; i++){
      fscanf(hyperfile, "%f \n", &hypers[i]);
      //printf("%f\n", hypers[i]);
      }
  model->x_mean=hypers[0];
  model->x_std=hypers[1];
  model->y_mean=hypers[2];
  model->y_std=hypers[3];
  model->grad_std=hypers[4];
  printf("xmean %f , xstd %f , ymean %f , ystd %f , grad_std %f \n", model->x_mean, model->x_std, model->y_mean, model->y_std, model->grad_std);
  }

void check_nn_status(tf_model *model){
  // make sure everything worked
  if(TF_GetCode(model->status) == TF_OK)
  {
    printf("TF status OK\n");
  }
  else
  {
    printf("%s", TF_Message(model->status));
  }
}


void input_data_to_model(tf_model *model, float data[][3], int natom){
  /* pass data to model */

  int64_t dims[] = {1,natom,DIM};
  int ndims = 3; //length of above array
  int ndata = sizeof(float)*natom*DIM; // This is tricky, it's number of bytes not number of elements

  /* To get a prediction now we allocate the input tensor with the input data */;

  TF_Tensor* int_tensor = TF_NewTensor(TF_FLOAT, dims, ndims, data, ndata, &NoOpDeallocator, 0);
  // make sure everything worked
  // if (int_tensor != NULL)
  //  {
  //     printf("TF_NewTensor is OK\n");
  //  }
  // else
  //  {
  //   printf("ERROR: Failed TF_NewTensor\n");
  //  }

  // and pass the tensor to the model
  model->input_val[0] = int_tensor;

}


void get_nn_hamiltonian(charge_transfer_t *ct, t_state *state_global, t_mdatoms *mdatoms, dftb_t *dftb)
{
/* NN Model trained with AIMAT group's code for learning energies & forces
 * Note: Model takes coordinates in Angstrom and returns energies in Hartree
 * and forces in Hartree/nm.
 * NOTE: We abuse the dftb struct's fields to pass the gradients to the rest of the code.
 * This avoids having to go through the entire codebase to implement ML force support.
 */
  long    iatom, isite, jsite, i, j, k, natoms;
  int atom_qm_idx, idx_i;
  long idx_lookup[ct->sites][ct->site[0].atoms];
  int oth_lookup[2*ct->site[0].atoms];
  float coords[ct->sites][ct->site[0].atoms][DIM];
  double tmp_frc, wf_factor;

  //  store the coordinates of atoms in QM zone
  //  and create lookup for indices needed for forces
  for (iatom = 0; iatom < mdatoms->homenr; iatom++)
  {
    atom_qm_idx = 0;
    for (i = 0; i < ct->sites; i++)
    {
      for (j = 0; j < ct->site[i].atoms; j++)
      {
        if (iatom == ct->site[i].atom[j])
        {
          for (k = 0; k < DIM; k++)
          {
            coords[i][j][k] = (float)state_global->x.rvec_array()[iatom][k]*10.0;
          }
          idx_lookup[i][j] = atom_qm_idx;
        }
        atom_qm_idx++;
      }
    }
  }

  // // faster version, but might be buggy for forces:
  // long iatom_global, site, iatom_local;
  // int m;
  //
  // atom_qm_idx = 0;
  // for (site = 0; site < ct->sites; site++)
  // {
  //   for (iatom_local = 0; iatom_local < ct->site[i].atoms; iatom_local++)
  //   {
  //     for (m = 0; m < DIM; m++)
  //     {
  //       iatom_global = ct->site[site].atom[iatom_local];
  //       coords[site][iatom_local][m] = (float)state_global->x[iatom_global][m]*10.0;
  //     }
  //     // construct lookup array for forces later
  //     idx_lookup[site][iatom_local] = atom_qm_idx;
  //     atom_qm_idx++;
  //   }
  // }

  natoms = 2*ct->site[0].atoms;

  //coordinate array for pairs
  float coordinates[natoms][DIM];

  // zero out all terms of the hamiltonian - why? MK
  for (isite = 0; isite < ct->sites; isite++)
      for (jsite = 0; jsite < ct->sites; jsite++)
          ct->hamiltonian_ml[isite][jsite] = 0.0;

  // off-diagonal terms first, because of scaling... MK
  natoms = 2*ct->site[0].atoms;
  // allocate array for forces (reused for every prediction)
  float** offdiag_forces;
  offdiag_forces = (float**) malloc(natoms*sizeof(float*));
  tf_model *coupling_model = ct->coupling_model;

  int pairs_done = 0;
  double offdiag_start =  (double)clock()/CLOCKS_PER_SEC;
  printf("calculating offdiag values start at %f\n", offdiag_start);

  double nl_start =  (double)clock()/CLOCKS_PER_SEC;

  if (ct->neighbors_only == 1 || (ct->neighbors_only==-1 && ct->first_step==1) || (ct->neighbors_only>1 && (ct->step % ct->neighbors_only)==0)){
    int counter1, counter2;
    double dist;
    dvec bond;
    // get neighbor list
    counter1=-1;
    for (int i=0; i<ct->sites; i++)
    for (int ii=0; ii<ct->site[i].atoms; ii++){
      counter1++;
      counter2=-1;
      for (int j=0; j<ct->sites; j++)
      for (int jj=0; jj<ct->site[j].atoms; jj++){
        counter2++;
        if (counter2 > counter1) {break;}
        dvec_sub(dftb->phase1[i].x[ii],dftb->phase1[j].x[jj], bond);
        dist = dnorm(bond);
        if (dist < 20.0*BOHR2NM){
          dftb->nl[counter1][counter2]=dftb->nl[counter2][counter1]=1;
        }else{
          dftb->nl[counter1][counter2]=dftb->nl[counter2][counter1]=0;
        }
      }
    }
  }

  double nl_end =  (double)clock()/CLOCKS_PER_SEC;
  printf("calculating neighborlist took %f\n", nl_end - nl_start);

  for (isite = 0; isite < ct->sites - 1; isite++)
  {
    for (jsite = isite + 1; jsite < ct->sites; jsite++)
    {
      if (ct->neighbors_only!=0 && (!dftb->nl[isite][jsite]))
      {
        continue;
      }

      pairs_done++;

      for (i = 0 ; i < ct->site[0].atoms; i++)
      {
        for (k = 0 ; k < DIM; k++)
        {
          coordinates[i][k] = (coords[isite][i][k] - coupling_model->x_mean)/coupling_model->x_std;
          coordinates[i+ct->site[0].atoms][k] = (coords[jsite][i][k] - coupling_model->x_mean)/coupling_model->x_std;
         // coordinates[i+ct->site[0].atoms][k] = coords[jsite][i][k]; no scaling
        }
        oth_lookup[i] = idx_lookup[isite][i];
        oth_lookup[i+ct->site[0].atoms] = idx_lookup[jsite][i];
      }

      input_data_to_model(coupling_model, coordinates, natoms);
      // Get the predictions
      TF_SessionRun(coupling_model->session, NULL, coupling_model->input, coupling_model->input_val, coupling_model->num_inputs, coupling_model->output, coupling_model->output_val, coupling_model->num_outputs, NULL, 0,NULL , coupling_model->status);
      //check_nn_status(coupling_model);

      float* buffer = (float*)TF_TensorData(coupling_model->output_val[0]);
      ct->hamiltonian_ml[isite][jsite] = buffer[0] * coupling_model->y_std + coupling_model->y_mean;
//          printf("coupling model predicted %f, that's %f Hartree\n", buffer[0], ct->hamiltonian_ml[isite][jsite]);
      ct->hamiltonian_ml[jsite][isite] = ct->hamiltonian_ml[isite][jsite];

      // we get the forces as 1D array (with natoms*3 elements)
      float* offsets_force = (float*) TF_TensorData(coupling_model->output_val[1]);
      // we reverse the Scaling

      for (int i=0; i<DIM*natoms; i++)
      {
        offsets_force[i] *= coupling_model->grad_std;
      }
      // we can reshape into 2D array in-place
      for(int i=0; i<natoms; i++)
      {
        offdiag_forces[i] = offsets_force + i*DIM;
      }

      // now we pass that to dftb->phase2.grad in the order we got earlier:
      for(int l=0; l < natoms; l++)
      {
        idx_i = oth_lookup[l];
        //printf("indices of atoms %d %d\n", idx_i, idx_j);
        for(int m=0; m<DIM; m++)
        {
          //MK tag forces: set wf_factor = 1 for pure gradients w/o wf
          if (ct->pure_forces){
            wf_factor = 1.0;
          }else if(ct->do_ncv){
            wf_factor = ct->tfs_vector[ct->ncv_isurf][isite]*ct->tfs_vector[ct->ncv_jsurf][jsite];
          }else{
            wf_factor = ct->wf[isite]*ct->wf[jsite] + ct->wf[isite+ct->dim]*ct->wf[jsite+ct->dim];
          }
          tmp_frc = wf_factor * offdiag_forces[l][m];
          dftb->phase2.grad[idx_i][m] += ct->is_hole_transfer ? tmp_frc : -1.0*tmp_frc;
        }
//        printf("frc for atom %d:  %f %f %f\n", idx_i, dftb->phase2.grad[idx_i][0], dftb->phase2.grad[idx_i][1], dftb->phase2.grad[idx_i][2]);
      }
    }
  }

  double offdiag_end = (double)clock()/CLOCKS_PER_SEC;
  double diff = offdiag_end - nl_end;

  printf("offdiagonal predictions took %f, for %d pairs %f per pair\n", diff, pairs_done, diff/(double)pairs_done);

  free(offdiag_forces);

  // set up site prediction
  natoms = ct->site[0].atoms;
  // allocate array for forces (reused for every prediction)
  float** diag_forces;
  diag_forces = (float**) malloc(natoms*sizeof(float*));

  tf_model *site_model = ct->site_model;

  //scale according to nn hyperparams
  for (i=0; i < ct->sites; i++)
    {
      for (j=0; j < ct->site[i].atoms; j++)
        {
          for (k=0; k<DIM; k++)
          {
            coords[i][j][k] = (coords[i][j][k] - site_model->x_mean)/site_model->x_std;
//              coords[i][j][k] = coords[i][j][k]; //no scaling
          }
        }
    }

  double diag_start = (double)clock()/CLOCKS_PER_SEC;
  printf("calculating diag values start at %f\n", diag_start);
  for (int isite = 0; isite < ct->sites; isite++)
  {
      input_data_to_model(site_model, coords[isite], natoms);
      // Get the predictions
      TF_SessionRun(site_model->session, NULL, site_model->input, site_model->input_val, site_model->num_inputs, site_model->output, site_model->output_val, site_model->num_outputs, NULL, 0,NULL , site_model->status);
      //check_nn_status(site_model);
      float* buffer = (float*)TF_TensorData(site_model->output_val[0]);
      ct->hamiltonian_ml[isite][isite] = buffer[0] * site_model->y_std + site_model->y_mean;
//      printf("site model predicted %f, that's %f Hartree\n", buffer[0], ct->hamiltonian_ml[isite][isite]);

      // we get the forces as 1D array (with natoms*3 elements)
      float* offsets_force = (float*) TF_TensorData(site_model->output_val[1]);
      // we reverse the Scaling
      for (int i=0; i<DIM*natoms; i++){
        offsets_force[i] = offsets_force[i] * site_model->grad_std;
      }
      // we can reshape into 2D array in-place
      for(int i=0; i<natoms; i++)
      {
        diag_forces[i] = offsets_force + i*DIM;
      }
      for(int i=0; i<natoms; i++)
      {
        for(int m=0; m<DIM; m++)
        {
          //MK tag forces: set wf_factor = 1 for pure gradients w/o wf
          if (ct->pure_forces){
            wf_factor = 1.0;
          }else if (ct->do_ncv){
            wf_factor = ct->tfs_vector[ct->ncv_isurf][ct->indFO[isite]]*ct->tfs_vector[ct->ncv_jsurf][ct->indFO[isite]];
          }else{
            wf_factor = (ct->wf[ct->indFO[isite]]*ct->wf[ct->indFO[isite]] + ct->wf[ct->indFO[isite]+ct->dim]*ct->wf[ct->indFO[isite]+ct->dim]);
          }
          dftb->phase1[isite].grad[i][m] = wf_factor*diag_forces[i][m];
        }
//        printf("diag_frc %f %f %f\n",  dftb->phase1[isite].grad[i][0],  dftb->phase1[isite].grad[i][1], dftb->phase1[isite].grad[i][2]);

      }
      if (ct->jobtype != cteNOMOVEMENT && ct->lambda != 0.0)
      {
      // implicit reorganization
      ct->occupation[isite] = SQR(ct->tfs_diab[isite]) + SQR(ct->tfs_diab[isite+ct->dim]);
      ct->hamiltonian_ml[isite][isite] -= ct->lambda*ct->occupation[isite];
      }
  }

  double diag_end = (double)clock()/CLOCKS_PER_SEC;
  printf("diagonal predictions took %f, %f per site\n", diag_end - diag_start,  (diag_end-diag_start)/(double)ct->sites);

  free(diag_forces);

  printf("hamiltonian calculation done\n");
}

#endif

