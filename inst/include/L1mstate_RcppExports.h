// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#ifndef RCPP_L1mstate_RCPPEXPORTS_H_GEN_
#define RCPP_L1mstate_RCPPEXPORTS_H_GEN_

#include <RcppEigen.h>
#include <Rcpp.h>

namespace L1mstate {

    using namespace Rcpp;

    namespace {
        void validateSignature(const char* sig) {
            Rcpp::Function require = Rcpp::Environment::base_env()["require"];
            require("L1mstate", Rcpp::Named("quietly") = true);
            typedef int(*Ptr_validate)(const char*);
            static Ptr_validate p_validate = (Ptr_validate)
                R_GetCCallable("L1mstate", "_L1mstate_RcppExport_validate");
            if (!p_validate(sig)) {
                throw Rcpp::function_not_exported(
                    "C++ function with signature '" + std::string(sig) + "' not found in L1mstate");
            }
        }
    }

    inline List scaleQ(Eigen::MatrixXd X) {
        typedef SEXP(*Ptr_scaleQ)(SEXP);
        static Ptr_scaleQ p_scaleQ = NULL;
        if (p_scaleQ == NULL) {
            validateSignature("List(*scaleQ)(Eigen::MatrixXd)");
            p_scaleQ = (Ptr_scaleQ)R_GetCCallable("L1mstate", "_L1mstate_scaleQ");
        }
        RObject rcpp_result_gen;
        {
            RNGScope RCPP_rngScope_gen;
            rcpp_result_gen = p_scaleQ(Shield<SEXP>(Rcpp::wrap(X)));
        }
        if (rcpp_result_gen.inherits("interrupted-error"))
            throw Rcpp::internal::InterruptedException();
        if (Rcpp::internal::isLongjumpSentinel(rcpp_result_gen))
            throw Rcpp::LongjumpException(rcpp_result_gen);
        if (rcpp_result_gen.inherits("try-error"))
            throw Rcpp::exception(Rcpp::as<std::string>(rcpp_result_gen).c_str());
        return Rcpp::as<List >(rcpp_result_gen);
    }

    inline double pletaQ(Eigen::VectorXd& xb, Eigen::VectorXi& nevent, Eigen::VectorXi& nevent1, Eigen::VectorXi& loc1, int& n) {
        typedef SEXP(*Ptr_pletaQ)(SEXP,SEXP,SEXP,SEXP,SEXP);
        static Ptr_pletaQ p_pletaQ = NULL;
        if (p_pletaQ == NULL) {
            validateSignature("double(*pletaQ)(Eigen::VectorXd&,Eigen::VectorXi&,Eigen::VectorXi&,Eigen::VectorXi&,int&)");
            p_pletaQ = (Ptr_pletaQ)R_GetCCallable("L1mstate", "_L1mstate_pletaQ");
        }
        RObject rcpp_result_gen;
        {
            RNGScope RCPP_rngScope_gen;
            rcpp_result_gen = p_pletaQ(Shield<SEXP>(Rcpp::wrap(xb)), Shield<SEXP>(Rcpp::wrap(nevent)), Shield<SEXP>(Rcpp::wrap(nevent1)), Shield<SEXP>(Rcpp::wrap(loc1)), Shield<SEXP>(Rcpp::wrap(n)));
        }
        if (rcpp_result_gen.inherits("interrupted-error"))
            throw Rcpp::internal::InterruptedException();
        if (Rcpp::internal::isLongjumpSentinel(rcpp_result_gen))
            throw Rcpp::LongjumpException(rcpp_result_gen);
        if (rcpp_result_gen.inherits("try-error"))
            throw Rcpp::exception(Rcpp::as<std::string>(rcpp_result_gen).c_str());
        return Rcpp::as<double >(rcpp_result_gen);
    }

    inline Eigen::VectorXd d1Q(Eigen::VectorXd& xb, Eigen::VectorXd& tevent, int& N, Eigen::VectorXi& nevent, Eigen::VectorXi& nevent1, Eigen::VectorXi& loc1, int& n) {
        typedef SEXP(*Ptr_d1Q)(SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP);
        static Ptr_d1Q p_d1Q = NULL;
        if (p_d1Q == NULL) {
            validateSignature("Eigen::VectorXd(*d1Q)(Eigen::VectorXd&,Eigen::VectorXd&,int&,Eigen::VectorXi&,Eigen::VectorXi&,Eigen::VectorXi&,int&)");
            p_d1Q = (Ptr_d1Q)R_GetCCallable("L1mstate", "_L1mstate_d1Q");
        }
        RObject rcpp_result_gen;
        {
            RNGScope RCPP_rngScope_gen;
            rcpp_result_gen = p_d1Q(Shield<SEXP>(Rcpp::wrap(xb)), Shield<SEXP>(Rcpp::wrap(tevent)), Shield<SEXP>(Rcpp::wrap(N)), Shield<SEXP>(Rcpp::wrap(nevent)), Shield<SEXP>(Rcpp::wrap(nevent1)), Shield<SEXP>(Rcpp::wrap(loc1)), Shield<SEXP>(Rcpp::wrap(n)));
        }
        if (rcpp_result_gen.inherits("interrupted-error"))
            throw Rcpp::internal::InterruptedException();
        if (Rcpp::internal::isLongjumpSentinel(rcpp_result_gen))
            throw Rcpp::LongjumpException(rcpp_result_gen);
        if (rcpp_result_gen.inherits("try-error"))
            throw Rcpp::exception(Rcpp::as<std::string>(rcpp_result_gen).c_str());
        return Rcpp::as<Eigen::VectorXd >(rcpp_result_gen);
    }

    inline Eigen::VectorXd d2Q(Eigen::VectorXd& xb, Eigen::VectorXd& tevent, int& N, Eigen::VectorXi& nevent, Eigen::VectorXi& nevent1, Eigen::VectorXi& loc1, int& n) {
        typedef SEXP(*Ptr_d2Q)(SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP);
        static Ptr_d2Q p_d2Q = NULL;
        if (p_d2Q == NULL) {
            validateSignature("Eigen::VectorXd(*d2Q)(Eigen::VectorXd&,Eigen::VectorXd&,int&,Eigen::VectorXi&,Eigen::VectorXi&,Eigen::VectorXi&,int&)");
            p_d2Q = (Ptr_d2Q)R_GetCCallable("L1mstate", "_L1mstate_d2Q");
        }
        RObject rcpp_result_gen;
        {
            RNGScope RCPP_rngScope_gen;
            rcpp_result_gen = p_d2Q(Shield<SEXP>(Rcpp::wrap(xb)), Shield<SEXP>(Rcpp::wrap(tevent)), Shield<SEXP>(Rcpp::wrap(N)), Shield<SEXP>(Rcpp::wrap(nevent)), Shield<SEXP>(Rcpp::wrap(nevent1)), Shield<SEXP>(Rcpp::wrap(loc1)), Shield<SEXP>(Rcpp::wrap(n)));
        }
        if (rcpp_result_gen.inherits("interrupted-error"))
            throw Rcpp::internal::InterruptedException();
        if (Rcpp::internal::isLongjumpSentinel(rcpp_result_gen))
            throw Rcpp::LongjumpException(rcpp_result_gen);
        if (rcpp_result_gen.inherits("try-error"))
            throw Rcpp::exception(Rcpp::as<std::string>(rcpp_result_gen).c_str());
        return Rcpp::as<Eigen::VectorXd >(rcpp_result_gen);
    }

    inline double max_lambdaQ(List& x) {
        typedef SEXP(*Ptr_max_lambdaQ)(SEXP);
        static Ptr_max_lambdaQ p_max_lambdaQ = NULL;
        if (p_max_lambdaQ == NULL) {
            validateSignature("double(*max_lambdaQ)(List&)");
            p_max_lambdaQ = (Ptr_max_lambdaQ)R_GetCCallable("L1mstate", "_L1mstate_max_lambdaQ");
        }
        RObject rcpp_result_gen;
        {
            RNGScope RCPP_rngScope_gen;
            rcpp_result_gen = p_max_lambdaQ(Shield<SEXP>(Rcpp::wrap(x)));
        }
        if (rcpp_result_gen.inherits("interrupted-error"))
            throw Rcpp::internal::InterruptedException();
        if (Rcpp::internal::isLongjumpSentinel(rcpp_result_gen))
            throw Rcpp::LongjumpException(rcpp_result_gen);
        if (rcpp_result_gen.inherits("try-error"))
            throw Rcpp::exception(Rcpp::as<std::string>(rcpp_result_gen).c_str());
        return Rcpp::as<double >(rcpp_result_gen);
    }

    inline List l1msQ(List& x, Eigen::VectorXd lambda, int nlambda, int p, double thresh, int maxit) {
        typedef SEXP(*Ptr_l1msQ)(SEXP,SEXP,SEXP,SEXP,SEXP,SEXP);
        static Ptr_l1msQ p_l1msQ = NULL;
        if (p_l1msQ == NULL) {
            validateSignature("List(*l1msQ)(List&,Eigen::VectorXd,int,int,double,int)");
            p_l1msQ = (Ptr_l1msQ)R_GetCCallable("L1mstate", "_L1mstate_l1msQ");
        }
        RObject rcpp_result_gen;
        {
            RNGScope RCPP_rngScope_gen;
            rcpp_result_gen = p_l1msQ(Shield<SEXP>(Rcpp::wrap(x)), Shield<SEXP>(Rcpp::wrap(lambda)), Shield<SEXP>(Rcpp::wrap(nlambda)), Shield<SEXP>(Rcpp::wrap(p)), Shield<SEXP>(Rcpp::wrap(thresh)), Shield<SEXP>(Rcpp::wrap(maxit)));
        }
        if (rcpp_result_gen.inherits("interrupted-error"))
            throw Rcpp::internal::InterruptedException();
        if (Rcpp::internal::isLongjumpSentinel(rcpp_result_gen))
            throw Rcpp::LongjumpException(rcpp_result_gen);
        if (rcpp_result_gen.inherits("try-error"))
            throw Rcpp::exception(Rcpp::as<std::string>(rcpp_result_gen).c_str());
        return Rcpp::as<List >(rcpp_result_gen);
    }

    inline List cvl1msQ(List& x, Eigen::VectorXd lambda, int nlambda, int p, double thresh, int maxit, List& xF) {
        typedef SEXP(*Ptr_cvl1msQ)(SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP);
        static Ptr_cvl1msQ p_cvl1msQ = NULL;
        if (p_cvl1msQ == NULL) {
            validateSignature("List(*cvl1msQ)(List&,Eigen::VectorXd,int,int,double,int,List&)");
            p_cvl1msQ = (Ptr_cvl1msQ)R_GetCCallable("L1mstate", "_L1mstate_cvl1msQ");
        }
        RObject rcpp_result_gen;
        {
            RNGScope RCPP_rngScope_gen;
            rcpp_result_gen = p_cvl1msQ(Shield<SEXP>(Rcpp::wrap(x)), Shield<SEXP>(Rcpp::wrap(lambda)), Shield<SEXP>(Rcpp::wrap(nlambda)), Shield<SEXP>(Rcpp::wrap(p)), Shield<SEXP>(Rcpp::wrap(thresh)), Shield<SEXP>(Rcpp::wrap(maxit)), Shield<SEXP>(Rcpp::wrap(xF)));
        }
        if (rcpp_result_gen.inherits("interrupted-error"))
            throw Rcpp::internal::InterruptedException();
        if (Rcpp::internal::isLongjumpSentinel(rcpp_result_gen))
            throw Rcpp::LongjumpException(rcpp_result_gen);
        if (rcpp_result_gen.inherits("try-error"))
            throw Rcpp::exception(Rcpp::as<std::string>(rcpp_result_gen).c_str());
        return Rcpp::as<List >(rcpp_result_gen);
    }

}

#endif // RCPP_L1mstate_RCPPEXPORTS_H_GEN_
