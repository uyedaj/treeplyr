#include <Rcpp.h>
#include <cmath>
#include <vector>

RcppExport SEXP geiger_descendants (SEXP phy)
{
    /* requires preorder (pruningwise ordering) of 'phylo' object */
  try {
  	Rcpp::List phylo(phy);
		int N = Rcpp::as<int>(phylo["N"]);
		int maxnode = Rcpp::as<int>(phylo["MAXNODE"]);
		std::vector<int> anc = Rcpp::as<std::vector<int> >(phylo["ANC"]);
		std::vector<int> des = Rcpp::as<std::vector<int> >(phylo["DES"]);
		
		int rows = maxnode-1;
		int root = N+1;
		
		std::vector< std::vector<int> > TIPS;
		std::vector< std::vector<int> > FDESC;
		std::vector< std::vector<int> > ADESC;
		std::vector< std::vector<int> > AANC;
        
		std::vector<int> empty;
		
		
		int i, j, k, s, t, z, dn, fd;
		
		/* initialize TIPS with known descendants (tips and root), otherwise leave empty */
		std::vector<int> cur;
		for(i = 0; i < maxnode; i++) {
			FDESC.push_back(empty);
            AANC.push_back(empty);
			if(i < N)
			{
				cur.push_back(i+1);
				TIPS.push_back(cur);
				ADESC.push_back(cur);
			}
			else
			{
				TIPS.push_back(empty);
				ADESC.push_back(empty);
			}
			cur.clear();
		}
		
		/* store nodes associated with root -- TIPS */
		for(i=0; i<N; i++){
			cur.push_back(i+1);
		}
		TIPS.at(N)=cur;
		
		/* store nodes associated with root -- ALL */
		cur.clear();
		for(i=0; i<maxnode; i++){
			if(i!=N)
			{
				cur.push_back(i+1);
			}
		}
		ADESC.at(N)=cur;
		
		/* store nodes associated with root -- FIRST */
		cur.clear();
		for(i=0; i<rows; i++){
			int idx = anc.at(i);
			if(idx==root)
			{
				cur.push_back(des.at(i));
			}
		}
		FDESC.at(N)=cur;
		
		/* collect descendants of each node in edge matrix (using pruningwise order to eliminate unnecessary computations) */
		for(i = 0; i < rows; i++){
			int nd = des.at(i);
			if(nd > N) {
				/* find descendant nodes */
				std::vector<int> subtends;
				for(j=0; j<rows; j++){
					int idx = anc.at(j);
					if(idx==nd)
					{
						subtends.push_back(des.at(j));
					}
				}
				FDESC.at(nd-1)=subtends;
				
				/* find immediate descendants of nd */
				std::vector<int> subtendedtips;
				std::vector<int> subtendednodes;
				s=subtends.size();
				for(k = 0; k < s; k++){
					
					/* find nodes subtended by immediate descendants of nd */
					fd = subtends.at(k);
					subtendednodes.push_back(fd);
                    
					if(fd<root) {
						subtendedtips.push_back(fd);
					} else {
						std::vector<int> descnodes = ADESC[fd-1];
						t=descnodes.size();
						for(z = 0; z < t; z++){
							dn=descnodes[z];
							if(dn<root)
							{
								subtendedtips.push_back(dn);
							}
							subtendednodes.push_back(dn);
						}
					}
				}
				/* store tips associated with nd into main list */
				TIPS.at(nd-1)=subtendedtips;
				ADESC.at(nd-1)=subtendednodes;
                s=subtendednodes.size();
                
                /* store ancestors */
                for(k=0; k<s; k++){
                    int idx = subtendednodes.at(k);
                    std::vector<int> ancnodes = AANC[idx-1];
                    ancnodes.push_back(nd);
                    AANC.at(idx-1)=ancnodes;
                }
			}
		}
		
		for(i=0; i<N; i++){
			ADESC.at(i)=empty;
		}
        for(k=0; k<maxnode; k++){
            int idx = k+1;
            if(idx!=root){
                std::vector<int> ancnodes = AANC[k];
                ancnodes.push_back(root);
                AANC.at(k)=ancnodes;
            }
        }
		return Rcpp::List::create(Rcpp::Named("tips",TIPS),
								  Rcpp::Named("fdesc",FDESC),
								  Rcpp::Named("adesc",ADESC),
                                  Rcpp::Named("anc",AANC));
	} catch( std::exception &ex ) {
		forward_exception_to_r( ex );
	} catch(...) {
		::Rf_error( "C++ exception: unknown reason" );
	}
	return R_NilValue;
}
