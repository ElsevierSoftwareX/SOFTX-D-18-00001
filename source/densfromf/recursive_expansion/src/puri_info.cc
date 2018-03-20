/* Ergo, version 3.6, a program for linear scaling electronic structure
 * calculations.
 * Copyright (C) 2017 Elias Rudberg, Emanuel H. Rubensson, Pawel Salek,
 * and Anastasia Kruchinina.
 * 
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 * 
 * Primary academic reference:
 * Kohn−Sham Density Functional Theory Electronic Structure Calculations 
 * with Linearly Scaling Computational Time and Memory Usage,
 * Elias Rudberg, Emanuel H. Rubensson, and Pawel Salek,
 * J. Chem. Theory Comput. 7, 340 (2011),
 * <http://dx.doi.org/10.1021/ct100611z>
 * 
 * For further information about Ergo, see <http://www.ergoscf.org>.
 */

/** @file puri_info.cc

    @brief File containing definitions of functions declared in classes IterationInfo and PuriInfo.
    IterationInfo is a class with the information stored for each iteration of the recursive expansion.
    PuriInfo is a class containing general information about the recursive expansion.

    @author Anastasia Kruchinina <em>responsible</em>
    @sa puri_info.h
*/



#include "puri_info.h"


typedef PuriInfo::real real;

void PuriInfo::get_poly_seq(std::vector<int> &poly)
{
  poly.clear();
  poly.resize(total_it+1);

  for(int it = 0; it <= total_it; ++it)
    poly[it] = Iterations[it].poly;
}

void PuriInfo::get_vec_frob_norms(std::vector<real> &norms)
{
  norms.clear();
  norms.resize(total_it+1);

  for(int it = 0; it <= total_it; ++it)
    norms[it] = Iterations[it].XmX2_fro_norm;
}


void PuriInfo::get_vec_mixed_norms(std::vector<real> &norms)
{
  norms.clear();
  norms.resize(total_it+1);

  for(int it = 0; it <= total_it; ++it)
    norms[it] = Iterations[it].XmX2_mixed_norm;
}


void PuriInfo::get_vec_traces(std::vector<real> & traces)
{
  traces.clear();
  traces.resize(total_it+1);

  for(int it = 0; it <= total_it; ++it)
    traces[it] = Iterations[it].XmX2_trace;
}


void PuriInfo::get_spectrum_bounds(real &lower_spectrum_bound_, real &upper_spectrum_bound_) const
{
  lower_spectrum_bound_ = lower_spectrum_bound;
  upper_spectrum_bound_ = upper_spectrum_bound;
}

void PuriInfo::set_spectrum_bounds(const real lower_spectrum_bound_, const real upper_spectrum_bound_)
{
  lower_spectrum_bound = lower_spectrum_bound_;
  upper_spectrum_bound = upper_spectrum_bound_;
}

void PuriInfo::print_collected_info()
{

  do_output(LOG_CAT_INFO, LOG_AREA_DENSFROMF,"=========================================================================================================================================");
  do_output(LOG_CAT_INFO, LOG_AREA_DENSFROMF,"                                           SHORT OUTPUT                                                ");
  do_output(LOG_CAT_INFO, LOG_AREA_DENSFROMF,"=========================================================================================================================================");

  // if used new stopping criterion
  if( stopping_criterion == 1)
    {
      do_output(LOG_CAT_INFO, LOG_AREA_DENSFROMF,"Used the new stopping criterion");
      do_output(LOG_CAT_INFO, LOG_AREA_DENSFROMF,"Estimated number of iterations:\t %d\n", estim_total_it);

      do_output(LOG_CAT_INFO, LOG_AREA_DENSFROMF,"Total number of iterations: iter + addit = %d + %d\n", total_it - additional_iterations, additional_iterations);
      do_output(LOG_CAT_INFO, LOG_AREA_DENSFROMF,"Iter \t Poly \t Order     Frob \t Trace \t\t NNZ(X)  \t NNZ(X2) \t  Treshold \t Total time");
      do_output(LOG_CAT_INFO, LOG_AREA_DENSFROMF,"------------------------------------------------------------------------------------------------------------------------");
      for(int i = 1; i <= total_it - additional_iterations; i++)
	{
	  do_output(LOG_CAT_INFO, LOG_AREA_DENSFROMF,"%3d \t %1d \t %6.2lf\t %12.5e \t %12.5e \t %6.2lf \t %6.2lf \t %9.5e \t %9.5lf", Iterations[i].it, Iterations[i].poly, Iterations[i].order,
	  	    (double)Iterations[i].XmX2_fro_norm, (double)Iterations[i].XmX2_trace, Iterations[i].NNZ_X, Iterations[i].NNZ_X2,
	  	    (double)Iterations[i].threshold_X, (double)Iterations[i].total_time);
	}
      do_output(LOG_CAT_INFO, LOG_AREA_DENSFROMF,"-----------------------------------------------------------------");
      for(int i = total_it - additional_iterations + 1; i <= total_it; i++)
	{
	  do_output(LOG_CAT_INFO, LOG_AREA_DENSFROMF,"%3d \t %1d \t %6.2lf\t %12.5e \t %12.5e \t %6.2lf \t %6.2lf \t %9.5e \t %9.5lf", Iterations[i].it, Iterations[i].poly, Iterations[i].order,
	  	    (double)Iterations[i].XmX2_fro_norm, (double)Iterations[i].XmX2_trace, Iterations[i].NNZ_X, Iterations[i].NNZ_X2,
	  	    (double)Iterations[i].threshold_X, (double)Iterations[i].total_time);
	}
      do_output(LOG_CAT_INFO, LOG_AREA_DENSFROMF,"Total time of the purification: \t %lf sec\n", (double)total_time);

      do_output(LOG_CAT_INFO, LOG_AREA_DENSFROMF,"");
      if(method == 2) // accelerated
	{
	  do_output(LOG_CAT_INFO, LOG_AREA_DENSFROMF,"Iter \t Poly \t  Alpha \t\t Lumo \t\t\t Homo \t\t C ");
	  do_output(LOG_CAT_INFO, LOG_AREA_DENSFROMF,"-----------------------------------------------------------------");
	  for(int i = 1; i <= total_it - additional_iterations; i++)
	    {
	      do_output(LOG_CAT_INFO, LOG_AREA_DENSFROMF,"%d \t %d \t %lf \t [%lf, %lf] \t [%lf, %lf] \t %lf \n", Iterations[i].it, Iterations[i].poly, (double)Iterations[i].alpha, (double)Iterations[i].lumo_bound_low, (double)Iterations[i].lumo_bound_upp, (double)Iterations[i].homo_bound_low, (double)Iterations[i].homo_bound_upp, (double)Iterations[i].constantC); 
	    }
	  do_output(LOG_CAT_INFO, LOG_AREA_DENSFROMF,"-----------------------------------------------------------------");
	  for(int i = total_it - additional_iterations + 1; i <= total_it; i++)
	    {
	      do_output(LOG_CAT_INFO, LOG_AREA_DENSFROMF,"%d \t %d \t %lf \t [%lf, %lf] \t [%lf, %lf] \t %lf \n", Iterations[i].it, Iterations[i].poly, (double)Iterations[i].alpha, (double)Iterations[i].lumo_bound_low, (double)Iterations[i].lumo_bound_upp, (double)Iterations[i].homo_bound_low, (double)Iterations[i].homo_bound_upp, (double)Iterations[i].constantC); 
	    }
	}
	  
	  
#ifdef CHECK_IF_STOPPED_TOO_LATE_OR_TOO_EARLY
      double XmX2_fro_norm_at_stop      = Iterations[total_it - additional_iterations    ].XmX2_fro_norm;
      double XmX2_fro_norm_before_stop  = Iterations[total_it - additional_iterations - 2].XmX2_fro_norm;
      double XmX2_fro_norm_after_stop   = Iterations[total_it - additional_iterations + 2].XmX2_fro_norm;
      if(XmX2_fro_norm_at_stop > XmX2_fro_norm_before_stop*2)
	do_output(LOG_CAT_INFO, LOG_AREA_DENSFROMF,"FAILED (STOPPED TOO LATE)! XmX2_fro_norm_before_stop = %e, XmX2_fro_norm_at_stop = %e.\n", XmX2_fro_norm_before_stop, XmX2_fro_norm_at_stop);
      if(XmX2_fro_norm_at_stop > XmX2_fro_norm_after_stop*2 && XmX2_fro_norm_after_stop != 0.0)
	do_output(LOG_CAT_INFO, LOG_AREA_DENSFROMF,"FAILED (STOPPED TOO EARLY)! XmX2_fro_norm_after_stop = %e, XmX2_fro_norm_at_stop = %e.\n", XmX2_fro_norm_after_stop, XmX2_fro_norm_at_stop);
#endif
    }
  else // old stopping criterion
    {
      do_output(LOG_CAT_INFO, LOG_AREA_DENSFROMF,"Used the old stopping criterion");
      do_output(LOG_CAT_INFO, LOG_AREA_DENSFROMF,"Estimated number of iterations:\t %d\n", estim_total_it);

      do_output(LOG_CAT_INFO, LOG_AREA_DENSFROMF,"Total number of iterations: iter + addit = %d + %d\n", total_it - additional_iterations, additional_iterations);
      do_output(LOG_CAT_INFO, LOG_AREA_DENSFROMF,"");
      do_output(LOG_CAT_INFO, LOG_AREA_DENSFROMF,"Iter \t Frob \t\t Trace \t\t  NNZ(X) NNZ(X2) Total time ");
      do_output(LOG_CAT_INFO, LOG_AREA_DENSFROMF,"-----------------------------------------------------------------------------");
      for(int i = 1; i <= total_it - additional_iterations; i++)
	{
	  do_output(LOG_CAT_INFO, LOG_AREA_DENSFROMF,"%d \t %e \t %e \t %d \t %d \t %lf \n", Iterations[i].it, 
		    (double)Iterations[i].XmX2_fro_norm, (double)Iterations[i].XmX2_trace, Iterations[i].NNZ_X, Iterations[i].NNZ_X2,
		    (double)Iterations[i].total_time); 
	}
      do_output(LOG_CAT_INFO, LOG_AREA_DENSFROMF,"");
      do_output(LOG_CAT_INFO, LOG_AREA_DENSFROMF,"Total time of the purification: \t %lf sec\n", (double)total_time);
    }
  do_output(LOG_CAT_INFO, LOG_AREA_DENSFROMF,"=========================================================================================================================================");

}



void PuriInfo::print_collected_info_printf()
{

  printf("=========================================================================================================================================\n");
  printf("=                                           SHORT OUTPUT                                                =\n");
  printf("=========================================================================================================================================\n");

  // if used new stopping criterion
  if( stopping_criterion == 1)
    {
      printf("Used the new stopping criterion\n");
      printf("Estimated number of iterations:\t %d\n\n", estim_total_it);

      printf("Total number of iterations: iter + addit = %d + %d\n", total_it - additional_iterations, additional_iterations);
      printf("Iter \t Poly \t Order     Frob \t Trace \t\t NNZ(X)  \t NNZ(X2) \t  Treshold \t Total time\n");
      printf("-------------------------------------------------------------------------------------------------------------------------------------\n");
      for(int i = 1; i <= total_it - additional_iterations; i++)
	{
	  printf("%3d \t %1d \t %6.2lf\t %12.5e \t %12.5e \t %6.2lf \t %6.2lf \t %9.5e \t %9.5lf\n", Iterations[i].it, 
		 Iterations[i].poly, (double)Iterations[i].order,
		 (double)Iterations[i].XmX2_fro_norm, (double)Iterations[i].XmX2_trace, (double)Iterations[i].NNZ_X, (double)Iterations[i].NNZ_X2,
		 (double)Iterations[i].threshold_X, (double)Iterations[i].total_time);
	}
      printf("-------------------------------------------------------------------------------------------------------------------------------------\n");
      for(int i = total_it - additional_iterations + 1; i <= total_it; i++)
	{
	  printf("%3d \t %1d \t %6.2lf\t %12.5e \t %12.5e \t %6.2lf \t %6.2lf \t %9.5e \t %9.5lf\n", Iterations[i].it, 
		 Iterations[i].poly, (double)Iterations[i].order,
		 (double)Iterations[i].XmX2_fro_norm, (double)Iterations[i].XmX2_trace, (double)Iterations[i].NNZ_X, (double)Iterations[i].NNZ_X2,
		 (double)Iterations[i].threshold_X, (double)Iterations[i].total_time);
	}
      printf("Total time of the purification: \t %lf sec\n", (double)total_time);

      printf("\n");
      if(method == 2) // accelerated
	{
	  printf("Iter \t Poly \t  Alpha \t\t Lumo \t\t\t Homo \t\t C \n");
	  printf("-----------------------------------------------------------------\n");
	  for(int i = 1; i <= total_it - additional_iterations; i++)
	    {
	      printf("%d \t %d \t %lf \t [%lf, %lf] \t [%lf, %lf] \t %lf \n", Iterations[i].it, Iterations[i].poly, 
		     (double)Iterations[i].alpha, (double)Iterations[i].lumo_bound_low, (double)Iterations[i].lumo_bound_upp, 
		     (double)Iterations[i].homo_bound_low, (double)Iterations[i].homo_bound_upp, (double)Iterations[i].constantC); 
	    }
	  printf("-------------------------------------------------------------------------------------------------------------------------------------\n");
	  for(int i = total_it - additional_iterations + 1; i <= total_it; i++)
	    {
	      printf("%d \t %d \t %lf \t [%lf, %lf] \t [%lf, %lf] \t %lf \n", Iterations[i].it, Iterations[i].poly, 
		     (double)Iterations[i].alpha, (double)Iterations[i].lumo_bound_low, (double)Iterations[i].lumo_bound_upp, 
		     (double)Iterations[i].homo_bound_low, (double)Iterations[i].homo_bound_upp, (double)Iterations[i].constantC); 
	    }
	}
	  
	  
#ifdef CHECK_IF_STOPPED_TOO_LATE_OR_TOO_EARLY
      double XmX2_fro_norm_at_stop      = Iterations[total_it - additional_iterations    ].XmX2_fro_norm;
      double XmX2_fro_norm_before_stop  = Iterations[total_it - additional_iterations - 2].XmX2_fro_norm;
      double XmX2_fro_norm_after_stop   = Iterations[total_it - additional_iterations + 2].XmX2_fro_norm;
      if(XmX2_fro_norm_at_stop > XmX2_fro_norm_before_stop*2)
	printf("FAILED (STOPPED TOO LATE)! XmX2_fro_norm_before_stop = %e, XmX2_fro_norm_at_stop = %e.\n", XmX2_fro_norm_before_stop, XmX2_fro_norm_at_stop);
      if(XmX2_fro_norm_at_stop > XmX2_fro_norm_after_stop*2 && XmX2_fro_norm_after_stop != 0.0)
	printf("FAILED (STOPPED TOO EARLY)! XmX2_fro_norm_after_stop = %e, XmX2_fro_norm_at_stop = %e.\n", XmX2_fro_norm_after_stop, XmX2_fro_norm_at_stop);
#endif
    }
  else // old stopping criterion
    {
      printf("Used the old stopping criterion\n");
      printf("Estimated number of iterations:\t %d\n", estim_total_it);

      printf("Total number of iterations: iter + addit = %d + %d\n", total_it - additional_iterations, additional_iterations);
      printf("\n");
      printf("Iter \t Frob \t\t Trace \t\t  NNZ(X) NNZ(X2) Total time \n");
      printf("-----------------------------------------------------------------------------\n");
      for(int i = 1; i <= total_it - additional_iterations; i++)
	{
	  printf("%d \t %e \t %e \t %.2lf \t %.2lf \t %lf \n", Iterations[i].it, 
		    (double)Iterations[i].XmX2_fro_norm, (double)Iterations[i].XmX2_trace, (double)Iterations[i].NNZ_X, (double)Iterations[i].NNZ_X2,
		    (double)Iterations[i].total_time); 
	}
      printf("\n");
      printf("Total time of the purification: \t %lf sec\n", (double)total_time);
    }
  printf("=========================================================================================================================================\n");

}

