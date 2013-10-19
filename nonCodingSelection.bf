function PopulateModelMatrix (dummy)
	{
for (h = 0; h<64; h=h+1) {if (_Genetic_Code[h]==10) {ModelMatrixDimension = ModelMatrixDimension-1;}} /* stop codon */
MG94HKY85_Matrix = {ModelMatrixDimension,ModelMatrixDimension};

hshift = 0;
for (h=0; h<64; h=h+1)
{
	if (_Genetic_Code[h]==10) 
	{
		hshift = hshift+1;
	}
	else
	{
		vshift = hshift;
		for (v = h+1; v<64; v=v+1)
		{
			diff = v-h;
			if (_Genetic_Code[v]==10) 
			{
				vshift = vshift+1;
			}
			else
			{
			  	if ((h$4==v$4)||((diff%4==0)&&(h$16==v$16))||(diff%16==0)) /* one step */
			  	{
			  		if (h$4==v$4)
			  		{
			  			transition = v%4;
			  			transition2= h%4;
			  		}
			  		else
			  		{
			  			if(diff%16==0)
			  			{
			  				transition = v$16;
			  				transition2= h$16;
			  			}
			  			else
			  			{
			  				transition = v%16$4;
			  				transition2= h%16$4;
			  			}
			  		}
			  		if (_Genetic_Code[0][h]==_Genetic_Code[0][v]) /* synonymous */
			  		{
			  			if (Abs(transition-transition2)%2) /* transversion */
			  			{
			  				MG94HKY85_Matrix[h-hshift][v-vshift] := kappa*t*ref_pis__[transition__];
			  				MG94HKY85_Matrix[v-vshift][h-hshift] := kappa*t*ref_pis__[transition2__];
			  			}
			  			else
			  			{
			  				MG94HKY85_Matrix[h-hshift][v-vshift] := t*ref_pis__[transition__];
			  				MG94HKY85_Matrix[v-vshift][h-hshift] := t*ref_pis__[transition2__];
			  			}
				  	}
			  		else
			  		{
			  			if (Abs(transition-transition2)%2) /* transversion */
			  			{
			  				MG94HKY85_Matrix[h-hshift][v-vshift] := omega*kappa*t*ref_pis__[transition__];
			  				MG94HKY85_Matrix[v-vshift][h-hshift] := omega*kappa*t*ref_pis__[transition2__];
			  			}
			  			else
			  			{
			  				MG94HKY85_Matrix[h-hshift][v-vshift] := omega*t*ref_pis__[transition__];
			  				MG94HKY85_Matrix[v-vshift][h-hshift] := omega*t*ref_pis__[transition2__];
			  			}
		  			}
			  	}
			 }
		 }
	}	
}

return 0;
	}

function BuildCodonFrequencies (obsF)
{
	PIStop = 1.0;
	result = {ModelMatrixDimension,1};
	hshift = 0;
	for (h=0; h<64; h=h+1)
		{
		first = h$16;
		second = h%16$4;
		third = h%4;
		if (_Genetic_Code[h]==10) 
			{
			hshift = hshift+1;
			/*PIStop = PIStop-obsF[first][0]*obsF[second][1]*obsF[third][2];*/
			PIStop = PIStop-ref_pis[first]*ref_pis[second]*ref_pis[third];
			continue; 
			}
		result[h-hshift]=ref_pis[first]*ref_pis[second]*ref_pis[third];
		/*result[h-hshift][0]=obsF[first][0]*obsF[second][1]*obsF[third][2];*/
	}
	return result*(1.0/PIStop);
}


function calculateEB(nsites,best_fs,d)
	{
	BEBFlag=0;
	NEBFlag=0;
	if ((EBKind==0) || (EBKind==2)) {BEBFlag=1;}
	if ((EBKind==1) || (EBKind==2)) {NEBFlag=1;}
	nClass=Columns(best_fs);
	if (NEBFlag==1)
		{
	  	NEB={nsites,nClass};
		fprintf(stdout,"calculate NEB\n");
		LikelihoodFunction query_L = (dsf_query, quer_hyphy_tree);
		ConstructCategoryMatrix(cm, query_L, COMPLETE); /*gets likelihood per site*/
		for (i=0; i<nsites;i=i+1)
			{
			marginal=0.0;
			for (k=0;k<nClass;k=k+1) {marginal=marginal+best_fs[1][k]*cm[k][i];}
			for (k=0;k<nClass;k=k+1) {NEB[i][k]=(best_fs[1][k]*cm[k][i])/marginal;}
			}
		}
	if (BEBFlag==1)
		{
		fprintf(stdout,"calculate BEB\n");
		fprintf(stdout,"discretization= ",d,"\n");
		BEB={nsites,nClass};
		/*d=number of discretization for the priors should be specified in the file that calls BEB*/
		ngrid=d*d*d*d; /*size of 4D grid*/
		fk={ngrid*nClass,1}; /*this will store f(zeta_k|n_s) for each point in 4D grid*/
		fhk={ngrid*nClass*npatt,1}; /*this will store f(xh|zeta(h)=zeta_k) for each point in 4D grid*/
		lnfXs={ngrid,1}; /*this log of the denominator*/
		postSite={npatt*nClass,1}; /*list of posterior probabilities for each site*/
		fns=(1/d)^4; /*prior f(ns)= flat priors in this case*/
		/*if (alternateKind==1) if this is M2A*/ 		/*To generalize according to the number of parameters!*/
		ZetaRange={2,2};
		ZetaRange[0][0]=0.0;
		ZetaRange[0][1]=0.95;
		ZetaRange[1][0]=1.05;
		ZetaRange[1][1]=11.0;
		GetDataInfo(duplicateMap,dsf_query); /*get the mapping of patterns to sites*/
		/*Get likelihood per site with parameters of the grid*/
		s=0; /*counter to store grid values in a vector*/
		for (y=0; y<d*d;y=y+1) /*loop through zeta fractions (=f0 f1 f2 0-1 in a ternary graph )*/
			{
			i=Sqrt(y)$1; /*i= f0 index of the triangle*/
			j=y-(i*i); /*j= f1 index of the triangle*/
			f0=(1.0+((j/2)$1)*3+(j % 2))/(3.0*d); /*f0 value from ternary graph as the center of the small triangle*/
			faux=(1.0+(d-1-i)*3+(j % 2))/(3.0*d);  /*f1 value from ternary graph as the center of the small triangle*/
			cc=0; /*counter for zeta2*/
			for (z=0; z<d*d;z=z+1) /*loop through zetas d*d modified*/
				{
				zeta0=(z$d)*((ZetaRange[0][1]-ZetaRange[0][0])/(d-1))+ZetaRange[0][0]; /*zeta0 0.0<=z0<=0.95*/
				zeta2=(cc)*((ZetaRange[1][1]-ZetaRange[1][0])/(d-1))+ZetaRange[1][0]; /*zeta2 1.05<=z0<=11.0*/
				cc=cc+1;
				if (cc>d-1) cc=0;
				/*here we keep branch length and kappa from previous optimization since they appear much less important to the
				calculation if the posterior probabilities, and their values are fixed at the MLEs (Yang et al.2005)*/
				LikelihoodFunction L = (dsf_query, quer_hyphy_tree);
				LFCompute (L,LF_START_COMPUTE); 
				LFCompute (L,res2);  /*recalculate likelihood with new f0, f1 f2 and zeta0, zeta2*/
				LFCompute (L,LF_DONE_COMPUTE);
				ConstructCategoryMatrix(cm, L, COMPLETE); /*gets likelihood per site*/
				GetInformation(distrInfo,site_class);
				/*this is an optimization. Instead of going through all sites, we are going through all site patterns, we will multiply
				by the number of sites per pattern at the end*/
				doneAlready={}; /*this will remember which sites have been "visited"*/
				for (h=0;h<nsites;h=h+1) /*goes through all sites*/
					{
					map2unique=duplicateMap[h]; /*get the pattern site index*/
					if (doneAlready[map2unique]==0)
						{
						doneAlready[map2unique]=1;
						for (k=0;k<nClass;k=k+1) {fhk[s*(nClass*npatt)+nClass*map2unique+k]=cm[k][h]; /*f(x_h|zeta(h)=zeta_k)*/}
						}
					}
				for (k=0;k<nClass;k=k+1) {fk[nClass*s+k]=distrInfo[1][k];} /*f(zeta_k|n_s)*/
				s=s+1;
				}
			}
			
		/* At this point we have calculated f(x_h|zeta(h)=zeta_k) for all site patterns, all zeta classes for all grid values*/
		/* We have also calculated f(zeta_k|n_s) for all grid values and all zeta classes*/
		/*Calculate Log(fx) by integrating over priors*/
		fx=0.0; /*denominator initialized, this will in fact be the log of f(x) -log are easier to handl eat thi spoint because of the small float problem*/
		for (s=0;s<ngrid;s=s+1) /*go through all grid points*/
			{
			lnfXs[s]=0.0;
			for (h=0;h<npatt;h=h+1) /*goes though all site patterns -here we going to used logs because of the following sumlogpexp to fix the small float problem*/
				{
				fh=0;
				for (k=0;k<nClass;k=k+1) {fh=fh+(fhk[s*(nClass*npatt)+nClass*h+k]*fk[nClass*s+k]);} /*gets f(xh) that is the sum of f(x_h|zeta(h)=zeta_k)*f(zeta_k|n_s) for all zeta k*/
				lnfXs[s]=lnfXs[s]+(Log(fh)*dsf_query.site_freqs[h]); /* denominator fx=product of f(xh) for all sites = sum of logs - here we multiply as well by the number of sites per partern upstream_data_set_filter.site_freqs[h]*/
				}
			lnfXs[s]=lnfXs[s]+Log(fns); /*we multiply by the prior =add log of the prior*/
			/*at this point we have the denominator f(x) for one point of the grid*/
			/*sumlogexp trick for small flaots while integrating fx for each grid point */
			if (s==0) /*first time, assign lnfXs[s] as the maximum value*/
				{
				scalingFactor=lnfXs[s];
				continue;
				}
			t=lnfXs[s]-scalingFactor;
			if (t>0) /*change scalingFactor because it is a new maximum*/
				{
				fx=fx*Exp(-t)+1;
				scalingFactor=lnfXs[s];
				}
			else
				{fx=fx+Exp(t);}
			}
		fx=Log(fx)+scalingFactor;;
		/*Integration of Log of Denominator f(x) completed*/
		
		/*Calculate posterior probability f(zeta|X)*/
		for (h=0;h<npatt;h=h+1) /*goes through all site patterns*/
			{
			for (thisClass=0;thisClass<nClass;thisClass=thisClass+1) /*go through all zeta classes*/
				{
				/*Now for a given site pattern and a given zeta class we are going to compute the BEB posterior probability of haveing this site pattern "h" belonging to the class "thisClass"*/
				for (s=0;s<ngrid;s=s+1) /*goes though all points of the grid for the integration*/
					{
					/*calculation of lnfx_h= ln(f(x_h|zeta(h)=zeta_k)*f(zeta_k|n_s)*(Product of fj for all sites j but h))
										   = ln(f(x_h|zeta(h)=zeta_k)*f(zeta_k|n_s))-ln(fh)+lnfXs[s]*/
					fh=0;
					for (k=0;k<nClass;k=k+1) {fh=fh+(fhk[s*(nClass*npatt)+nClass*h+k]*fk[nClass*s+k]);}
					fhSite=(fhk[s*(nClass*npatt)+nClass*h+thisClass]*fk[nClass*s+thisClass])/fh;
					lnfxs_h=Log(fhSite)+lnfXs[s];
					/*once again we are using the sumlogexp trick and that's why we used logs just above*/
					if ((thisClass==0) && (s==0)) /*first time, assign lnfxs_h as the maximum value*/
						{
						scalingFactor=lnfxs_h;
						continue;
						}
					if (lnfxs_h>scalingFactor) /*change scalingFactor because it is a new maximum*/
						{
						for (k=0;k<nClass;k=k+1) {postSite[k*npatt+h]=postSite[k*npatt+h]*Exp(scalingFactor-lnfxs_h);}
						scalingFactor=lnfxs_h;
						}
					postSite[thisClass*npatt+h]=postSite[thisClass*npatt+h]+Exp(lnfxs_h-scalingFactor); /*this is the integration of the numerator*/
					}
				}
			for (k=0;k<nClass;k=k+1) {postSite[k*npatt+h]=postSite[k*npatt+h]*Exp(scalingFactor-fx);} /*numerator/denominator= numerator*Exp(-Log(f(x)) after numerator being scaled back to numerator*exp(scalingfactor) for the sumlogexp trick*/
			}
		/*Now we are done!*/
		map2unique=duplicateMap[h];
		for (h=0;h<nsites;h=h+1) 
			{
			map2unique=duplicateMap[h];
			for (k=0;k<nClass;k=k+1)
				{
				BEB[h][k]=postSite[k*npatt+map2unique];
				}
			}
		}
	return 0;
	}

function setNeutral(dummy)
	{
	if (neutralProxy==1)  /*use HKY85*/
		{
		ref_mat = {{*, kappa*t, t, kappa*t}
           		{kappa*t, *, kappa*t, t}
           		{t, kappa*t, *, kappa*t}
           		{kappa*t, t, kappa*t, *}};
		Model ref_submod = (ref_mat, ref_pis, 1);
		}
	else  /*use MG94xHKY85*/
		{
     	global omega = 1;
     	ModelMatrixDimension = 64;
    	PopulateModelMatrix(0);

    	codonFreqs = BuildCodonFrequencies (ref_pis);
  		/*define the codon model */
		Model ref_submod = (MG94HKY85_Matrix,codonFreqs,0);
		}
	return 0;
	}

function fitNullModel(dummy)
	{
	/* Set up query compartment submodel. and reference compartment submodel.*/
	global kappa;
	quer_mat = {{*, kappa*t, t, kappa*t}
				{kappa*t, *, kappa*t, t}
				{t, kappa*t, *, kappa*t}
				{kappa*t, t, kappa*t, *}};
	Model quer_submod = (quer_mat, quer_pis, 1);
	setNeutral(0);
	if (nullKind==0) /*null1*/
		{
		log_Ls = {fit_repl_count, 1};
		for (repl_i = 0; repl_i < fit_repl_count; repl_i = repl_i+1)
			{
			kappa = Random(0.0, 10.0);
			UseModel(quer_submod);
			Tree quer_hyphy_tree = treeString;
			UseModel(ref_submod);
			Tree ref_hyphy_tree = treeString;
			ReplicateConstraint("this1.?.t := this2.?.t", quer_hyphy_tree, ref_hyphy_tree);
			LikelihoodFunction L = (dsf_ref, ref_hyphy_tree,dsf_query, quer_hyphy_tree);
			Optimize(MLEs, L);
			log_L = MLEs[1][0];
			log_Ls[repl_i] = log_L;
			if (repl_i == 0 || best_log_L < log_L)
				{
				best_repl_i = repl_i;
				best_log_L = log_L;
				best_MLEs = MLEs;
				}
			}
		estimated_kappa = best_MLEs[0][0];
		fprintf(stdout,"*** Null model ***","\n");
		fprintf(stdout,"inverse kappa: ",estimated_kappa,"\n");
		}
	else /*null2*/
		{
		global f0;
		f0 :< 1;
		global f_aux;
		f_aux :< 1;
		site_class_freqs = {{f0, (1.0-f0)*f_aux, (1.0-f0)*(1.0-f_aux)}};
		site_class_vals = {{0, 1, 2}};
		category site_class = (3, site_class_freqs, MEAN, , site_class_vals, 0, 2);
		global zeta0;
		zeta0 :< 1;
		global zeta_bgrnd;
		zeta_bgrnd := ((site_class == 0)+(site_class == 2))*zeta0 + (site_class == 1);
		global zeta_fgrnd;
		zeta_fgrnd := (site_class == 0)*zeta0 + ((site_class == 1)+(site_class == 2));
		log_Ls = {fit_repl_count, 1};
		for (repl_i = 0; repl_i < fit_repl_count; repl_i = repl_i+1)
			{
  			kappa = Random(0.0, 10.0);
  			f0 = Random(0.0, 1.0);
  			f_aux = Random(0.0, 1.0);
  			zeta0 = Random(0.0, 1.0);
  			UseModel(quer_submod);
  			Tree quer_hyphy_tree = treeString;
  			UseModel(ref_submod);
  			Tree ref_hyphy_tree = treeString;
  			ReplicateConstraint("this1.?.t := zeta_bgrnd*this2.?.t", quer_hyphy_tree, ref_hyphy_tree);
  			ExecuteCommands("quer_hyphy_tree."+fgrnd_branch_name+".t := zeta_fgrnd*ref_hyphy_tree."+fgrnd_branch_name+".t;");
  			LikelihoodFunction L = (dsf_query, quer_hyphy_tree, dsf_ref, ref_hyphy_tree);
  			Optimize(MLEs, L);
  			log_L = MLEs[1][0];
  			log_Ls[repl_i] = log_L;
  			if (repl_i == 0 || best_log_L < log_L)
  				{
    			best_repl_i = repl_i;
    			best_log_L = log_L;
    			best_MLEs = MLEs;
  				}
			}
		estimated_kappa = best_MLEs[0][3];
		estimated_f0 = best_MLEs[0][0];
		estimated_f_aux = best_MLEs[0][1];
		estimated_f1 = (1.0-estimated_f0)*estimated_f_aux;
		estimated_f2 = (1.0-estimated_f0)*(1.0-estimated_f_aux);
		estimated_zeta0 = best_MLEs[0][2];
		fprintf(stdout,"*** Null model ***","\n");
		fprintf(stdout,"inverse kappa: ",estimated_kappa,"\n");
		fprintf(stdout, "f0: ", estimated_f0,"\n");
		fprintf(stdout, "f1: ", estimated_f1,"\n");
		fprintf(stdout, "f2: ", estimated_f2,"\n");
		fprintf(stdout, "zeta0: ", estimated_zeta0,"\n");
		fprintf(stdout, "zeta1: 1","\n\n");
		}
	return best_MLEs;
	}

function fitAlternateModel(dummy)
	{
	/* Set up query compartment submodel and reference compartment submodel. */
	global kappa;	
	quer_mat = {{*, kappa*t, t, kappa*t}
            	{kappa*t, *, kappa*t, t}
            	{t, kappa*t, *, kappa*t}
            	{kappa*t, t, kappa*t, *}};
	Model quer_submod = (quer_mat, quer_pis, 1);
	setNeutral(0);
	if (alternateKind==0) /*alt1 no branch*/
		{
		global zeta;
		zeta :> 1;
		log_Ls = {fit_repl_count, 1};
		for (repl_i = 0; repl_i < fit_repl_count; repl_i = repl_i+1)
				{
  				kappa = Random(0.0, 10.0);
  				zeta = Random(1.0, 10.0);
  				UseModel(quer_submod);
  				Tree quer_hyphy_tree = treeString;
  				UseModel(ref_submod);
  				Tree ref_hyphy_tree = treeString;
  				ReplicateConstraint("this1.?.t := zeta*this2.?.t", quer_hyphy_tree, ref_hyphy_tree);
  				LikelihoodFunction L = (dsf_query, quer_hyphy_tree, dsf_ref, ref_hyphy_tree);
  				Optimize(MLEs, L);
  				log_L = MLEs[1][0];
  				log_Ls[repl_i] = log_L;
  				if (repl_i == 0 || best_log_L < log_L)
  						{
    					best_repl_i = repl_i;
    					best_log_L = log_L;
    					best_MLEs = MLEs;
  						}
  				}
  		estimated_kappa = best_MLEs[0][1];
		estimated_zeta = best_MLEs[0][0];
		fprintf(stdout,"*** Alternate model ***","\n");
		fprintf(stdout, "inverse kappa: ", estimated_kappa,"\n");
		fprintf(stdout, "zeta: ", estimated_zeta,"\n\n");
		}
	else
		{
		if (alternateKind==1) /*alt2 no branch*/
			{
			global f0;
			f0 :< 1;
			global f_aux;
			f_aux :< 1;
			site_class_freqs = {{f0, (1.0-f0)*f_aux, (1.0-f0)*(1.0-f_aux)}};
			site_class_vals = {{0, 1, 2}};
			category site_class = (3, site_class_freqs, MEAN, , site_class_vals, 0, 2);
			global zeta0;
			zeta0 :< 1;
			global zeta2;
			zeta2 :> 1;
			global zeta;
			zeta := (site_class == 0)*zeta0 + (site_class == 1) + (site_class == 2)*zeta2;
			log_Ls = {fit_repl_count, 1};
			for (repl_i = 0; repl_i < fit_repl_count; repl_i = repl_i+1)
				{
				kappa = Random(0.0, 10.0);
				f0 = Random(0.0, 1.0);
				f_aux = Random(0.0, 1.0);
				zeta0 = Random(0.0, 1.0);
				zeta2 = Random(1.0, 10.0);
				UseModel(quer_submod);
				Tree quer_hyphy_tree = treeString;
				UseModel(ref_submod);
				Tree ref_hyphy_tree = treeString;
				ReplicateConstraint("this1.?.t := zeta*this2.?.t", quer_hyphy_tree, ref_hyphy_tree);
				LikelihoodFunction L = (dsf_query, quer_hyphy_tree, dsf_ref, ref_hyphy_tree);
				Optimize(MLEs, L);
				log_L = MLEs[1][0];
				log_Ls[repl_i] = log_L;
				if (repl_i == 0 || best_log_L < log_L)
					{
					best_repl_i = repl_i;
					best_log_L = log_L;
					best_MLEs = MLEs;
					}
				}
		fprintf(stdout,"*** Alternate model ***","\n");
		estimated_kappa = best_MLEs[0][4];
		estimated_f0 = best_MLEs[0][0];
		estimated_f_aux = best_MLEs[0][1];
		estimated_f1 = (1.0-estimated_f0)*estimated_f_aux;
		estimated_f2 = (1.0-estimated_f0)*(1.0-estimated_f_aux);
		estimated_f3 = (1.0-estimated_f0)*(1.0-estimated_f_aux)*(1.0-estimated_f0)*estimated_f_aux/(estimated_f0 + (1.0-estimated_f0)*estimated_f_aux);
		estimated_zeta0 = best_MLEs[0][2];
		estimated_zeta2 = best_MLEs[0][3];
		fprintf(stdout, "f0: ", estimated_f0,"\n");
		fprintf(stdout, "f1: ", estimated_f1,"\n");
		fprintf(stdout, "f2: ", estimated_f2,"\n");
		fprintf(stdout, "f3: ", estimated_f3,"\n");
		fprintf(stdout, "zeta0: ", estimated_zeta0,"\n");
		fprintf(stdout, "zeta1: 1","\n");
		fprintf(stdout, "zeta2: ", estimated_zeta2,"\n\n");
		}
		else
			{
			if (alternateKind==2) /*alt1  branch spec*/
				{
				global zeta_fgrnd;
				zeta_fgrnd :> 1;
				/* Fit model to data sets. */
				log_Ls = {fit_repl_count, 1};
				for (repl_i = 0; repl_i < fit_repl_count; repl_i = repl_i+1)
					{
  					kappa = Random(0.0, 10.0);
  					zeta_fgrnd = Random(1.0, 10.0);
  					UseModel(quer_submod);
  					Tree quer_hyphy_tree = treeString;
  					UseModel(ref_submod);
  					Tree ref_hyphy_tree = treeString;
  					ReplicateConstraint("this1.?.t := this2.?.t", quer_hyphy_tree, ref_hyphy_tree);
  					ExecuteCommands("quer_hyphy_tree."+fgrnd_branch_name+".t := zeta_fgrnd*ref_hyphy_tree."+fgrnd_branch_name+".t;");
  					LikelihoodFunction L = (dsf_query, quer_hyphy_tree, dsf_ref, ref_hyphy_tree);
  					Optimize(MLEs, L);
  					log_L = MLEs[1][0];
  					log_Ls[repl_i] = log_L;
  					if (repl_i == 0 || best_log_L < log_L)
  						{
    					best_repl_i = repl_i;
    					best_log_L = log_L;
    					best_MLEs = MLEs;
  						}
					}
				fprintf(stdout,"*** Alternate model ***","\n");
				estimated_kappa = best_MLEs[0][1];
				estimated_zeta = best_MLEs[0][0];
				fprintf(stdout,"inverse kappa: ",estimated_kappa,"\n");
				fprintf(stdout, "zeta: ", estimated_zeta,"\n\n");
				}
			else /*alt2 branch spec*/
				{
				global f0;
				f0 :< 1;
				global f_aux;
				f_aux :< 1;
				site_class_freqs = {{f0, (1.0-f0)*f_aux, (1.0-f0)*(1.0-f_aux)*f0/(f0 + (1.0-f0)*f_aux), (1.0-f0)*(1.0-f_aux)*(1.0-f0)*f_aux/(f0 + (1.0-f0)*f_aux)}};
				site_class_vals = {{0, 1, 2, 3}};
				category site_class = (4, site_class_freqs, MEAN, , site_class_vals, 0, 3);
				global zeta0;
				zeta0 :< 1;
				global zeta2;
				zeta2 :> 1;
				global zeta_bgrnd;
				zeta_bgrnd := ((site_class == 0)+(site_class == 2))*zeta0 + ((site_class == 1)+(site_class == 3));
				global zeta_fgrnd;
				zeta_fgrnd := (site_class == 0)*zeta0 + (site_class == 1) + ((site_class == 2)+(site_class == 3))*zeta2;
				log_Ls = {fit_repl_count, 1};
				for (repl_i = 0; repl_i < fit_repl_count; repl_i = repl_i+1)
					{
				 	kappa = Random(0.0, 10.0);
				  	f0 = Random(0.0, 1.0);
				  	f_aux = Random(0.0, 1.0);
				  	zeta0 = Random(0.0, 1.0);
				  	zeta2 = Random(1.0, 10.0);
				  	UseModel(quer_submod);
				  	Tree quer_hyphy_tree = treeString;
				  	UseModel(ref_submod);
				  	Tree ref_hyphy_tree = treeString;
				  	ReplicateConstraint("this1.?.t := zeta_bgrnd*this2.?.t", quer_hyphy_tree, ref_hyphy_tree);
				  	ExecuteCommands("quer_hyphy_tree."+fgrnd_branch_name+".t := zeta_fgrnd*ref_hyphy_tree."+fgrnd_branch_name+".t;");
				 	LikelihoodFunction L = (dsf_query, quer_hyphy_tree, dsf_ref, ref_hyphy_tree);
				  	Optimize(MLEs, L);
				  	log_L = MLEs[1][0];
				  	log_Ls[repl_i] = log_L;
				  	if (repl_i == 0 || best_log_L < log_L)
				  		{
						best_repl_i = repl_i;
						best_log_L = log_L;
						best_MLEs = MLEs;
				  		}
					}
				/*alt2 branch spec*/
				estimated_kappa = best_MLEs[0][4];
				estimated_f0 = best_MLEs[0][0];
				estimated_f_aux = best_MLEs[0][1];
				estimated_f1 = (1.0-estimated_f0)*estimated_f_aux;
				estimated_f2 = (1.0-estimated_f0)*(1.0-estimated_f_aux);
				estimated_f3 = (1.0-estimated_f0)*(1.0-estimated_f_aux)*(1.0-estimated_f0)*estimated_f_aux/(estimated_f0 + (1.0-estimated_f0)*estimated_f_aux);
				estimated_zeta0 = best_MLEs[0][2];
				estimated_zeta2 = best_MLEs[0][3];
				fprintf(stdout,"*** Alternate model ***","\n");
				fprintf(stdout,"inverse kappa: ", estimated_kappa,"\n");
				fprintf(stdout, "f0: ", estimated_f0,"\n");
				fprintf(stdout, "f1: ", estimated_f1,"\n");
				fprintf(stdout, "f2: ", estimated_f2,"\n");
				fprintf(stdout, "f3: ", estimated_f3,"\n");
				fprintf(stdout, "zeta0: ", estimated_zeta0,"\n");
				fprintf(stdout, "zeta1: 1","\n");
				fprintf(stdout, "zeta2: ", estimated_zeta2,"\n\n");
				}
			}
		}
	return best_MLEs;
	}
	

function printSite(vectorClass)
	{
	maxClass=0;
	maxClassIndex=0;
	for (k=0;k<3;k=k+1)
		{
		if (vectorClass[k]>=maxClass)
			{
			maxClass=vectorClass[k];
			maxClassIndex=k;
			}
		}
	for (k=0;k<3;k=k+1)
		{
		fprintf(stdout,Format(vectorClass[k],10,6));
		if (k==maxClassIndex) {fprintf(stdout,"*");} else {fprintf(stdout," ");}
		}
	fprintf(stdout,"|");
	return 0;
	}

function printInformation(neutralProxy,branchSp,fgrnd_branch_name,nullKind,alternateKind,fit_repl_count)
	{
	fprintf (stdout,"\n");
	fprintf (stdout,"number of replicates for likelihood fitting= ",fit_repl_count,"\n");
	if (neutralProxy==0) {fprintf (stdout,"Coding region used as neutral proxy with MG94HKY85 model.\n");}
	if (neutralProxy==1) {fprintf (stdout,"Intron(s) used as neutral proxy with HKY85 model.\n");}
	if (branchSp==0)
		{
		fprintf (stdout,"Non-branch specific test.\n");
		fprintf (stdout,"NULL=Negative selection and Neutral evolution.\n");
		if (alternateKind==0) {fprintf (stdout,"ALTERNATE= positive selection.\n");}
		if (alternateKind==1) {fprintf (stdout,"ALTERNATE= Class 1: Negative, Class 2: neutral evolution and Class 3: positive selection.\n");}
		}
	if (branchSp==1)
		{
		fprintf (stdout,"Branch specific test with foreground branch ",fgrnd_branch_name,"\n");
		if (nullkind==0) {fprintf (stdout,"NULL= Class 1: Negative selection in FG and BG. Class 2: Neutral evolution in FG and BG.\n");}
		if (nullkind==1) {fprintf (stdout,"NULL= Class 1: Negative selection in FG and BG. Class 2: Neutral evolution in FG and BG. Class 3: Negative selection in BG, Neutral in FG.\n");}
		if (alternateKind==2) {fprintf (stdout,"ALTERNATE= Class 1: Negative selection or neutral in FG and BG. Class 2: Neutral or neutral evolution in BG and Positive in FG.\n");}
		if (alternateKind==3) {fprintf (stdout,"ALTERNATE= Class 1: Negative selection in FG and BG. Class 2: Neutral evolution in FG and BG. Class 3: Negative selection in BG, Positive in FG. Class 4: Neutral evolution in BG, Positive in FG.\n");}
		}
	fprintf (stdout,"\n");
	return 0;
	}
			
function printEB(dummy)
	{
	BEBFlag=0;
	NEBFlag=0;
	if ((EBKind==0) || (EBKind==2)) {BEBFlag=1;}
	if ((EBKind==1) || (EBKind==2)) {NEBFlag=1;}
	
	fprintf (stdout,"\n\n|--------|");
	if (BEBFlag==1) {fprintf (stdout,"---------------------------------|");}
	if (NEBFlag==1) {fprintf (stdout,"---------------------------------|");}
	fprintf (stdout,"\n|  Site  |");
	if (BEBFlag==1) {fprintf (stdout,"               BEB               |");}
	if (NEBFlag==1) {fprintf (stdout,"               NEB               |");}
	fprintf (stdout,"\n|        |");
	if (BEBFlag==1) {fprintf (stdout,"    Neg      Neutral      Pos    |");}
	if (NEBFlag==1) {fprintf (stdout,"    Neg      Neutral      Pos    |");}
	fprintf (stdout,"\n|--------|");
	if (BEBFlag==1) {fprintf (stdout,"---------------------------------|");}
	if (NEBFlag==1) {fprintf (stdout,"---------------------------------|");}
	

	for (h=0;h<nsites;h=h+1)
		{
		thisSite=h+1;
		fprintf(stdout,"\n| ",Format(thisSite,4,0),"   |");
		if (BEBFlag==1)
			{
			if (nClass==3) {printSite({{BEB[h][0],BEB[h][1],BEB[h][2]}});} else {printSite({{BEB[h][0],BEB[h][1],BEB[h][2]+BEB[h][3]}});}
			}
		if (NEBFlag==1)
			{
			if (nClass==3) {printSite({{NEB[h][0],NEB[h][1],NEB[h][2]}});} else {printSite({{NEB[h][0],NEB[h][1],NEB[h][2]+NEB[h][3]}});}
			}
		}
	fprintf (stdout,"\n|--------|");
	if (BEBFlag==1) {fprintf (stdout,"---------------------------------|");}
	if (NEBFlag==1) {fprintf (stdout,"---------------------------------|");}
	fprintf (stdout,"\n");
	return 0;
    }
/*--------------------------------------------------------------------------------------------------------------------------------------------------------------------*/
NICETY_LEVEL = 3;
LIKELIHOOD_FUNCTION_OUTPUT = 5;
MAXIMUM_ITERATIONS_PER_VARIABLE = 10000;
fit_repl_count=10;
d=10; /*BEB discretization*/

fprintf (stdout,"\nA tool for detecting positive selection on promoter sequences\n");
fprintf (stdout,"Implementation by Olivier Fedrigo (ofedrigo@duke.edu) and Ralph Haygood (rhaygood@duke,edu) 2007\n");
fprintf (stdout,"---------------------------------\n");



/*GET THE NEUTRAL PROXY*/
ChoiceList (neutralProxy, "Neutral proxy", 1, SKIP_NONE,"Synonymous sites", "Use synonymous substitution from a coding alignment","Intronic sites", "Use intronic sequences",);
if (neutralProxy < 0) {return 0;}
if (neutralProxy==0)
     {
     incFileName = HYPHY_LIB_DIRECTORY+"TemplateBatchFiles"+DIRECTORY_SEPARATOR+"TemplateModels"+DIRECTORY_SEPARATOR+"chooseGeneticCode.def";
     ExecuteCommands  ("#include \""+incFileName+"\";");
     SetDialogPrompt ("Please locate a coding alignment:");
     DataSet 	  ds_ref	   = ReadDataFile (PROMPT_FOR_FILE);
     DataSetFilter dsf_ref   = CreateFilter (ds_ref,3,"","",GeneticCodeExclusions);
     }
else
     {
     SetDialogPrompt ("Please locate an intronic alignment:");
     DataSet 	  ds_ref	   = ReadDataFile (PROMPT_FOR_FILE);
     DataSetFilter dsf_ref   = CreateFilter (ds_ref,1);
     }
neutral_path = LAST_FILE_PATH;
fprintf (stdout, "\nLoaded a ", dsf_ref.species, " sequence alignment with ", dsf_ref.sites, " nucleotides from\n",neutral_path,"\n");
HarvestFrequencies(ref_pis, ds_ref, 1, 1, 1);

/*GET THE TREE*/
incFileName = HYPHY_LIB_DIRECTORY+"TemplateBatchFiles"+DIRECTORY_SEPARATOR+"queryTree.bf";
ExecuteCommands  ("#include \""+incFileName+"\";");
ChoiceList (branchSp, "Branch Specific?", 1, SKIP_NONE,"No", "Wong and Nielsen test", "Yes", "Branch specific Wong and Nielsen test",);
if (branchSp < 0) {return 0;}
if (branchSp==1) /*branch specific*/
    {
	Tree givenTree = treeString;

    leafNodes	  = TipCount   (givenTree);	
    internalNodes = BranchCount(givenTree);
    choiceMatrix = {internalNodes+leafNodes,2};
    for (bc=0; bc<internalNodes; bc=bc+1)
        {
        choiceMatrix[bc][0] = BranchName(givenTree,bc);
        choiceMatrix[bc][1] = "Internal Branch Rooting " + givenTree[bc];
        }
    for (bc=0; bc<leafNodes; bc=bc+1)
        {
        choiceMatrix[bc+internalNodes][0] = TipName(givenTree,bc);
        choiceMatrix[bc+internalNodes][1] = "Leaf node " + choiceMatrix[bc+internalNodes][0];
        }
    ChoiceList (stOption,"Choose the foreground branch",1,SKIP_NONE,choiceMatrix);
    if (stOption < 0) {return 0;}
    fgrnd_branch_name= choiceMatrix[stOption][0];
    }

/*GET THE QUERY REGION*/
SetDialogPrompt ("Please locate a query alignment:");
DataSet ds_query  = ReadDataFile (PROMPT_FOR_FILE);
DataSetFilter dsf_query = CreateFilter (ds_query,1);
query_path = LAST_FILE_PATH;
fprintf (stdout, "\nLoaded a ", dsf_query.species, " sequence alignment with ", dsf_query.sites, " nucleotides from\n",query_path,"\n");
nsites=dsf_query.sites; /*number of sites*/
npatt=dsf_query.unique_sites; /*number of unique patterns*/
HarvestFrequencies(quer_pis, ds_query, 1, 1, 1);

/*CHECK IF DATASET ARE COMPATIBLE*/
if (ds_query.species != ds_ref.species)
    {
	fprintf (stdout, "\nError: both data sets must contain the same number of sequences\n");
	return 0;
    }
for (sid=0; sid<ds_coding.species; sid=sid+1)
    {
	GetString (neutralName, ds_ref, sid);
	GetString (queryName, ds_query, sid);	
	if (neutralName!=queryName)
	    {
		fprintf (stdout, "\nError: sequence name mismatch in sequence ", sid+1,"\nHad ", neutralName, " and ", queryName, "\n");
		return 0;
	    }
    }

/*GET MODEL COMPARISON TYPE*/
if (branchSp==0) /*no branch specific*/
    {
    nullKind=0;
    ChoiceList (alternateKind, "Alternate Model", 1, SKIP_NONE,"Null-Alternate 1", "NULL=Negative selection and Neutral evolution. ALTERNATE= positive selection","Null-Alternate 2", "NULL=Negative selection and Neutral evolution. ALTERNATE= Class 1: Negative, Class 2: neutral evolution and Class 3: positive selection.",);
    if (alternateKind < 0) {return 0;}

    }
else /*branch specific*/
    {
    ChoiceList (alternateKind, "Model comparisons", 1, SKIP_NONE,
					   "Null1-Alternate1", "NULL= Class 1: Negative selection in FG and BG. Class 2: Neutral evolution in FG and BG. ALTERNATE= Class 1: Negative selection or neutral in FG and BG. Class 2: Neutral or neutral evolution in BG and Positive in FG",
					   "Null2-Alternate2", "NULL= Class 1: Negative selection in FG and BG. Class 2: Neutral evolution in FG and BG. Class 3: Negative selection in BG, Neutral in FG. ALTERNATE= Class 1: Negative selection in FG and BG. Class 2: Neutral evolution in FG and BG. Class 3: Negative selection in BG, Positive in FG. Class 4: Neutral evolution in BG, Positive in FG.",);				   
    if (alternateKind < 0) {return 0;}
    nullKind=alternateKind;
    alternateKind=alternateKind+2;
	}

fprintf(stdout,alternateKind);

printInformation(neutralProxy,branchSp,fgrnd_branch_name,nullKind,alternateKind,fit_repl_count);

/*FIT MODELS*/
res_null=fitNullModel(0);
res_alt=fitAlternateModel(0);
fprintf(stdout,"Lk null = ",res_null[1][0]," Lk alt = ",res_alt[1][0],"\n");
LRT= -2*(res_null[1][0]-res_alt[1][0]);
fprintf(stdout,"LRT = ",LRT,"\n");
degFDiff = res_alt[1][1]-res_null[1][1];
pvalue = 1-CChi2(LRT,degFDiff);

fprintf(stdout,"LRT p-value = ",Format(pvalue,10,10),"\n");

if ((alternateKind==1) || (alternateKind==3))
{
ChoiceList (EBKind, "Empirical Bayes", 1, SKIP_NONE, "BEB", "Bayes Empirical Bayes","NEB", "Naive Empirical Bayes","Both", "Bayes Empirical Bayes and Naive Empirical Bayes",);

if (EBKind < 0) {return 0;}
GetInformation(best_fs,site_class); /*gets fraction for each zeta*/
nClass=Columns(best_fs);
NEB={nsites,nClass};
BEB={nsites,nClass};
calculateEB(nsites,best_fs,d);
if (BEBFlag==1)
	{
	fullSitesBEB    = {nsites,3};
	for (h=0;h<nsites;h=h+1)
		{
		if (nClass==3) 
			{
			neg=BEB[h][0];
			neutr=BEB[h][1];
			pos=BEB[h][2];
			fullSitesBEB[h][0]=neg;
			fullSitesBEB[h][1]=neutr;
			fullSitesBEB[h][2]=pos;
			}
		else
			{
			neg=BEB[h][0];
			neutr=BEB[h][1];
			pos=BEB[h][2]+BEB[h][3];
			fullSitesBEB[h][0]=neg;
			fullSitesBEB[h][1]=neutr;
			fullSitesBEB[h][2]=pos;
			}
		}
	labels = {{"Negative selection","Neutral evolution","Positive selection"}};
	/*Creates a windows with the output*/
	OpenWindow (CHARTWINDOW,{{"BEB"}
							   {"labels"},
							   {"fullSitesBEB"},
							   {"Bar Chart"},
							   {"Index"},
							   {labels[2]},
							   {"Site Index"},
							   {""},
							   {labels[2]},
							   {"0"}},
							   "SCREEN_WIDTH-60;SCREEN_HEIGHT-50;30;50");
	}
if (NEBFlag==1)
	{
	fullSitesNEB    = {nsites,3};
	for (h=0;h<nsites;h=h+1)
		{
		if (nClass==3) 
			{
			neg=NEB[h][0];
			neutr=NEB[h][1];
			pos=NEB[h][2];
			fullSitesNEB[h][0]=neg;
			fullSitesNEB[h][1]=neutr;
			fullSitesNEB[h][2]=pos;
			}
		else
			{
			neg=NEB[h][0];
			neutr=NEB[h][1];
			pos=NEB[h][2]+NEB[h][3];
			fullSitesNEB[h][0]=neg;
			fullSitesNEB[h][1]=neutr;
			fullSitesNEB[h][2]=pos;
			}
		}
	
	labels = {{"Negative selection","Neutral evolution","Positive selection"}};
	/*Creates a windows with the output*/
	OpenWindow (CHARTWINDOW,{{"NEB"}
							   {"labels"},
							   {"fullSitesNEB"},
							   {"Bar Chart"},
							   {"Index"},
							   {labels[2]},
							   {"Site Index"},
							   {""},
							   {labels[2]},
							   {"0"}},
							   "SCREEN_WIDTH-60;SCREEN_HEIGHT-50;30;50");
	}
printEB(0);
}

