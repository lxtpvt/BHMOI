// Generated by rstantools.  Do not edit by hand.

/*
    docbiBHM is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    docbiBHM is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with docbiBHM.  If not, see <http://www.gnu.org/licenses/>.
*/
#ifndef MODELS_HPP
#define MODELS_HPP
#define STAN__SERVICES__COMMAND_HPP
#include <rstan/rstaninc.hpp>
// Code generated by Stan version 2.21.0
#include <stan/model/model_header.hpp>
namespace model_tran_bhm_general_namespace {
using std::istream;
using std::string;
using std::stringstream;
using std::vector;
using stan::io::dump;
using stan::math::lgamma;
using stan::model::prob_grad;
using namespace stan::math;
static int current_statement_begin__;
stan::io::program_reader prog_reader__() {
    stan::io::program_reader reader;
    reader.add_event(0, 0, "start", "model_tran_bhm_general");
    reader.add_event(33, 31, "end", "model_tran_bhm_general");
    return reader;
}
#include <stan_meta_header.hpp>
class model_tran_bhm_general
  : public stan::model::model_base_crtp<model_tran_bhm_general> {
private:
        int K;
        std::vector<int> Y;
        std::vector<int> n;
        double mu0;
        double tau0;
        double a;
        double b;
public:
    model_tran_bhm_general(stan::io::var_context& context__,
        std::ostream* pstream__ = 0)
        : model_base_crtp(0) {
        ctor_body(context__, 0, pstream__);
    }
    model_tran_bhm_general(stan::io::var_context& context__,
        unsigned int random_seed__,
        std::ostream* pstream__ = 0)
        : model_base_crtp(0) {
        ctor_body(context__, random_seed__, pstream__);
    }
    void ctor_body(stan::io::var_context& context__,
                   unsigned int random_seed__,
                   std::ostream* pstream__) {
        typedef double local_scalar_t__;
        boost::ecuyer1988 base_rng__ =
          stan::services::util::create_rng(random_seed__, 0);
        (void) base_rng__;  // suppress unused var warning
        current_statement_begin__ = -1;
        static const char* function__ = "model_tran_bhm_general_namespace::model_tran_bhm_general";
        (void) function__;  // dummy to suppress unused var warning
        size_t pos__;
        (void) pos__;  // dummy to suppress unused var warning
        std::vector<int> vals_i__;
        std::vector<double> vals_r__;
        local_scalar_t__ DUMMY_VAR__(std::numeric_limits<double>::quiet_NaN());
        (void) DUMMY_VAR__;  // suppress unused var warning
        try {
            // initialize data block variables from context__
            current_statement_begin__ = 2;
            context__.validate_dims("data initialization", "K", "int", context__.to_vec());
            K = int(0);
            vals_i__ = context__.vals_i("K");
            pos__ = 0;
            K = vals_i__[pos__++];
            check_greater_or_equal(function__, "K", K, 0);
            current_statement_begin__ = 3;
            validate_non_negative_index("Y", "K", K);
            context__.validate_dims("data initialization", "Y", "int", context__.to_vec(K));
            Y = std::vector<int>(K, int(0));
            vals_i__ = context__.vals_i("Y");
            pos__ = 0;
            size_t Y_k_0_max__ = K;
            for (size_t k_0__ = 0; k_0__ < Y_k_0_max__; ++k_0__) {
                Y[k_0__] = vals_i__[pos__++];
            }
            current_statement_begin__ = 4;
            validate_non_negative_index("n", "K", K);
            context__.validate_dims("data initialization", "n", "int", context__.to_vec(K));
            n = std::vector<int>(K, int(0));
            vals_i__ = context__.vals_i("n");
            pos__ = 0;
            size_t n_k_0_max__ = K;
            for (size_t k_0__ = 0; k_0__ < n_k_0_max__; ++k_0__) {
                n[k_0__] = vals_i__[pos__++];
            }
            current_statement_begin__ = 5;
            context__.validate_dims("data initialization", "mu0", "double", context__.to_vec());
            mu0 = double(0);
            vals_r__ = context__.vals_r("mu0");
            pos__ = 0;
            mu0 = vals_r__[pos__++];
            current_statement_begin__ = 6;
            context__.validate_dims("data initialization", "tau0", "double", context__.to_vec());
            tau0 = double(0);
            vals_r__ = context__.vals_r("tau0");
            pos__ = 0;
            tau0 = vals_r__[pos__++];
            check_greater_or_equal(function__, "tau0", tau0, 0);
            current_statement_begin__ = 7;
            context__.validate_dims("data initialization", "a", "double", context__.to_vec());
            a = double(0);
            vals_r__ = context__.vals_r("a");
            pos__ = 0;
            a = vals_r__[pos__++];
            check_greater_or_equal(function__, "a", a, 0);
            current_statement_begin__ = 8;
            context__.validate_dims("data initialization", "b", "double", context__.to_vec());
            b = double(0);
            vals_r__ = context__.vals_r("b");
            pos__ = 0;
            b = vals_r__[pos__++];
            check_greater_or_equal(function__, "b", b, 0);
            // initialize transformed data variables
            // execute transformed data statements
            // validate transformed data
            // validate, set parameter ranges
            num_params_r__ = 0U;
            param_ranges_i__.clear();
            current_statement_begin__ = 12;
            validate_non_negative_index("p_k", "K", K);
            num_params_r__ += K;
            current_statement_begin__ = 13;
            num_params_r__ += 1;
            current_statement_begin__ = 14;
            num_params_r__ += 1;
        } catch (const std::exception& e) {
            stan::lang::rethrow_located(e, current_statement_begin__, prog_reader__());
            // Next line prevents compiler griping about no return
            throw std::runtime_error("*** IF YOU SEE THIS, PLEASE REPORT A BUG ***");
        }
    }
    ~model_tran_bhm_general() { }
    void transform_inits(const stan::io::var_context& context__,
                         std::vector<int>& params_i__,
                         std::vector<double>& params_r__,
                         std::ostream* pstream__) const {
        typedef double local_scalar_t__;
        stan::io::writer<double> writer__(params_r__, params_i__);
        size_t pos__;
        (void) pos__; // dummy call to supress warning
        std::vector<double> vals_r__;
        std::vector<int> vals_i__;
        current_statement_begin__ = 12;
        if (!(context__.contains_r("p_k")))
            stan::lang::rethrow_located(std::runtime_error(std::string("Variable p_k missing")), current_statement_begin__, prog_reader__());
        vals_r__ = context__.vals_r("p_k");
        pos__ = 0U;
        validate_non_negative_index("p_k", "K", K);
        context__.validate_dims("parameter initialization", "p_k", "vector_d", context__.to_vec(K));
        Eigen::Matrix<double, Eigen::Dynamic, 1> p_k(K);
        size_t p_k_j_1_max__ = K;
        for (size_t j_1__ = 0; j_1__ < p_k_j_1_max__; ++j_1__) {
            p_k(j_1__) = vals_r__[pos__++];
        }
        try {
            writer__.vector_lub_unconstrain(0, 1, p_k);
        } catch (const std::exception& e) {
            stan::lang::rethrow_located(std::runtime_error(std::string("Error transforming variable p_k: ") + e.what()), current_statement_begin__, prog_reader__());
        }
        current_statement_begin__ = 13;
        if (!(context__.contains_r("mu")))
            stan::lang::rethrow_located(std::runtime_error(std::string("Variable mu missing")), current_statement_begin__, prog_reader__());
        vals_r__ = context__.vals_r("mu");
        pos__ = 0U;
        context__.validate_dims("parameter initialization", "mu", "double", context__.to_vec());
        double mu(0);
        mu = vals_r__[pos__++];
        try {
            writer__.scalar_unconstrain(mu);
        } catch (const std::exception& e) {
            stan::lang::rethrow_located(std::runtime_error(std::string("Error transforming variable mu: ") + e.what()), current_statement_begin__, prog_reader__());
        }
        current_statement_begin__ = 14;
        if (!(context__.contains_r("tau")))
            stan::lang::rethrow_located(std::runtime_error(std::string("Variable tau missing")), current_statement_begin__, prog_reader__());
        vals_r__ = context__.vals_r("tau");
        pos__ = 0U;
        context__.validate_dims("parameter initialization", "tau", "double", context__.to_vec());
        double tau(0);
        tau = vals_r__[pos__++];
        try {
            writer__.scalar_lb_unconstrain(0, tau);
        } catch (const std::exception& e) {
            stan::lang::rethrow_located(std::runtime_error(std::string("Error transforming variable tau: ") + e.what()), current_statement_begin__, prog_reader__());
        }
        params_r__ = writer__.data_r();
        params_i__ = writer__.data_i();
    }
    void transform_inits(const stan::io::var_context& context,
                         Eigen::Matrix<double, Eigen::Dynamic, 1>& params_r,
                         std::ostream* pstream__) const {
      std::vector<double> params_r_vec;
      std::vector<int> params_i_vec;
      transform_inits(context, params_i_vec, params_r_vec, pstream__);
      params_r.resize(params_r_vec.size());
      for (int i = 0; i < params_r.size(); ++i)
        params_r(i) = params_r_vec[i];
    }
    template <bool propto__, bool jacobian__, typename T__>
    T__ log_prob(std::vector<T__>& params_r__,
                 std::vector<int>& params_i__,
                 std::ostream* pstream__ = 0) const {
        typedef T__ local_scalar_t__;
        local_scalar_t__ DUMMY_VAR__(std::numeric_limits<double>::quiet_NaN());
        (void) DUMMY_VAR__;  // dummy to suppress unused var warning
        T__ lp__(0.0);
        stan::math::accumulator<T__> lp_accum__;
        try {
            stan::io::reader<local_scalar_t__> in__(params_r__, params_i__);
            // model parameters
            current_statement_begin__ = 12;
            Eigen::Matrix<local_scalar_t__, Eigen::Dynamic, 1> p_k;
            (void) p_k;  // dummy to suppress unused var warning
            if (jacobian__)
                p_k = in__.vector_lub_constrain(0, 1, K, lp__);
            else
                p_k = in__.vector_lub_constrain(0, 1, K);
            current_statement_begin__ = 13;
            local_scalar_t__ mu;
            (void) mu;  // dummy to suppress unused var warning
            if (jacobian__)
                mu = in__.scalar_constrain(lp__);
            else
                mu = in__.scalar_constrain();
            current_statement_begin__ = 14;
            local_scalar_t__ tau;
            (void) tau;  // dummy to suppress unused var warning
            if (jacobian__)
                tau = in__.scalar_lb_constrain(0, lp__);
            else
                tau = in__.scalar_lb_constrain(0);
            // transformed parameters
            current_statement_begin__ = 18;
            validate_non_negative_index("theta", "K", K);
            Eigen::Matrix<local_scalar_t__, Eigen::Dynamic, 1> theta(K);
            stan::math::initialize(theta, DUMMY_VAR__);
            stan::math::fill(theta, DUMMY_VAR__);
            // transformed parameters block statements
            current_statement_begin__ = 19;
            stan::math::assign(theta, logit(p_k));
            // validate transformed parameters
            const char* function__ = "validate transformed params";
            (void) function__;  // dummy to suppress unused var warning
            current_statement_begin__ = 18;
            size_t theta_j_1_max__ = K;
            for (size_t j_1__ = 0; j_1__ < theta_j_1_max__; ++j_1__) {
                if (stan::math::is_uninitialized(theta(j_1__))) {
                    std::stringstream msg__;
                    msg__ << "Undefined transformed parameter: theta" << "(" << j_1__ << ")";
                    stan::lang::rethrow_located(std::runtime_error(std::string("Error initializing variable theta: ") + msg__.str()), current_statement_begin__, prog_reader__());
                }
            }
            // model body
            current_statement_begin__ = 24;
            lp_accum__.add(normal_log<propto__>(theta, mu, pow(tau, -(0.5))));
            current_statement_begin__ = 25;
            lp_accum__.add(normal_log<propto__>(mu, mu0, pow(tau0, -(0.5))));
            current_statement_begin__ = 26;
            lp_accum__.add(gamma_log<propto__>(tau, a, b));
            current_statement_begin__ = 28;
            for (int j = 1; j <= K; ++j) {
                current_statement_begin__ = 29;
                lp_accum__.add(binomial_log<propto__>(get_base1(Y, j, "Y", 1), get_base1(n, j, "n", 1), get_base1(p_k, j, "p_k", 1)));
            }
        } catch (const std::exception& e) {
            stan::lang::rethrow_located(e, current_statement_begin__, prog_reader__());
            // Next line prevents compiler griping about no return
            throw std::runtime_error("*** IF YOU SEE THIS, PLEASE REPORT A BUG ***");
        }
        lp_accum__.add(lp__);
        return lp_accum__.sum();
    } // log_prob()
    template <bool propto, bool jacobian, typename T_>
    T_ log_prob(Eigen::Matrix<T_,Eigen::Dynamic,1>& params_r,
               std::ostream* pstream = 0) const {
      std::vector<T_> vec_params_r;
      vec_params_r.reserve(params_r.size());
      for (int i = 0; i < params_r.size(); ++i)
        vec_params_r.push_back(params_r(i));
      std::vector<int> vec_params_i;
      return log_prob<propto,jacobian,T_>(vec_params_r, vec_params_i, pstream);
    }
    void get_param_names(std::vector<std::string>& names__) const {
        names__.resize(0);
        names__.push_back("p_k");
        names__.push_back("mu");
        names__.push_back("tau");
        names__.push_back("theta");
    }
    void get_dims(std::vector<std::vector<size_t> >& dimss__) const {
        dimss__.resize(0);
        std::vector<size_t> dims__;
        dims__.resize(0);
        dims__.push_back(K);
        dimss__.push_back(dims__);
        dims__.resize(0);
        dimss__.push_back(dims__);
        dims__.resize(0);
        dimss__.push_back(dims__);
        dims__.resize(0);
        dims__.push_back(K);
        dimss__.push_back(dims__);
    }
    template <typename RNG>
    void write_array(RNG& base_rng__,
                     std::vector<double>& params_r__,
                     std::vector<int>& params_i__,
                     std::vector<double>& vars__,
                     bool include_tparams__ = true,
                     bool include_gqs__ = true,
                     std::ostream* pstream__ = 0) const {
        typedef double local_scalar_t__;
        vars__.resize(0);
        stan::io::reader<local_scalar_t__> in__(params_r__, params_i__);
        static const char* function__ = "model_tran_bhm_general_namespace::write_array";
        (void) function__;  // dummy to suppress unused var warning
        // read-transform, write parameters
        Eigen::Matrix<double, Eigen::Dynamic, 1> p_k = in__.vector_lub_constrain(0, 1, K);
        size_t p_k_j_1_max__ = K;
        for (size_t j_1__ = 0; j_1__ < p_k_j_1_max__; ++j_1__) {
            vars__.push_back(p_k(j_1__));
        }
        double mu = in__.scalar_constrain();
        vars__.push_back(mu);
        double tau = in__.scalar_lb_constrain(0);
        vars__.push_back(tau);
        double lp__ = 0.0;
        (void) lp__;  // dummy to suppress unused var warning
        stan::math::accumulator<double> lp_accum__;
        local_scalar_t__ DUMMY_VAR__(std::numeric_limits<double>::quiet_NaN());
        (void) DUMMY_VAR__;  // suppress unused var warning
        if (!include_tparams__ && !include_gqs__) return;
        try {
            // declare and define transformed parameters
            current_statement_begin__ = 18;
            validate_non_negative_index("theta", "K", K);
            Eigen::Matrix<double, Eigen::Dynamic, 1> theta(K);
            stan::math::initialize(theta, DUMMY_VAR__);
            stan::math::fill(theta, DUMMY_VAR__);
            // do transformed parameters statements
            current_statement_begin__ = 19;
            stan::math::assign(theta, logit(p_k));
            if (!include_gqs__ && !include_tparams__) return;
            // validate transformed parameters
            const char* function__ = "validate transformed params";
            (void) function__;  // dummy to suppress unused var warning
            // write transformed parameters
            if (include_tparams__) {
                size_t theta_j_1_max__ = K;
                for (size_t j_1__ = 0; j_1__ < theta_j_1_max__; ++j_1__) {
                    vars__.push_back(theta(j_1__));
                }
            }
            if (!include_gqs__) return;
        } catch (const std::exception& e) {
            stan::lang::rethrow_located(e, current_statement_begin__, prog_reader__());
            // Next line prevents compiler griping about no return
            throw std::runtime_error("*** IF YOU SEE THIS, PLEASE REPORT A BUG ***");
        }
    }
    template <typename RNG>
    void write_array(RNG& base_rng,
                     Eigen::Matrix<double,Eigen::Dynamic,1>& params_r,
                     Eigen::Matrix<double,Eigen::Dynamic,1>& vars,
                     bool include_tparams = true,
                     bool include_gqs = true,
                     std::ostream* pstream = 0) const {
      std::vector<double> params_r_vec(params_r.size());
      for (int i = 0; i < params_r.size(); ++i)
        params_r_vec[i] = params_r(i);
      std::vector<double> vars_vec;
      std::vector<int> params_i_vec;
      write_array(base_rng, params_r_vec, params_i_vec, vars_vec, include_tparams, include_gqs, pstream);
      vars.resize(vars_vec.size());
      for (int i = 0; i < vars.size(); ++i)
        vars(i) = vars_vec[i];
    }
    std::string model_name() const {
        return "model_tran_bhm_general";
    }
    void constrained_param_names(std::vector<std::string>& param_names__,
                                 bool include_tparams__ = true,
                                 bool include_gqs__ = true) const {
        std::stringstream param_name_stream__;
        size_t p_k_j_1_max__ = K;
        for (size_t j_1__ = 0; j_1__ < p_k_j_1_max__; ++j_1__) {
            param_name_stream__.str(std::string());
            param_name_stream__ << "p_k" << '.' << j_1__ + 1;
            param_names__.push_back(param_name_stream__.str());
        }
        param_name_stream__.str(std::string());
        param_name_stream__ << "mu";
        param_names__.push_back(param_name_stream__.str());
        param_name_stream__.str(std::string());
        param_name_stream__ << "tau";
        param_names__.push_back(param_name_stream__.str());
        if (!include_gqs__ && !include_tparams__) return;
        if (include_tparams__) {
            size_t theta_j_1_max__ = K;
            for (size_t j_1__ = 0; j_1__ < theta_j_1_max__; ++j_1__) {
                param_name_stream__.str(std::string());
                param_name_stream__ << "theta" << '.' << j_1__ + 1;
                param_names__.push_back(param_name_stream__.str());
            }
        }
        if (!include_gqs__) return;
    }
    void unconstrained_param_names(std::vector<std::string>& param_names__,
                                   bool include_tparams__ = true,
                                   bool include_gqs__ = true) const {
        std::stringstream param_name_stream__;
        size_t p_k_j_1_max__ = K;
        for (size_t j_1__ = 0; j_1__ < p_k_j_1_max__; ++j_1__) {
            param_name_stream__.str(std::string());
            param_name_stream__ << "p_k" << '.' << j_1__ + 1;
            param_names__.push_back(param_name_stream__.str());
        }
        param_name_stream__.str(std::string());
        param_name_stream__ << "mu";
        param_names__.push_back(param_name_stream__.str());
        param_name_stream__.str(std::string());
        param_name_stream__ << "tau";
        param_names__.push_back(param_name_stream__.str());
        if (!include_gqs__ && !include_tparams__) return;
        if (include_tparams__) {
            size_t theta_j_1_max__ = K;
            for (size_t j_1__ = 0; j_1__ < theta_j_1_max__; ++j_1__) {
                param_name_stream__.str(std::string());
                param_name_stream__ << "theta" << '.' << j_1__ + 1;
                param_names__.push_back(param_name_stream__.str());
            }
        }
        if (!include_gqs__) return;
    }
}; // model
}  // namespace
typedef model_tran_bhm_general_namespace::model_tran_bhm_general stan_model;
#ifndef USING_R
stan::model::model_base& new_model(
        stan::io::var_context& data_context,
        unsigned int seed,
        std::ostream* msg_stream) {
  stan_model* m = new stan_model(data_context, seed, msg_stream);
  return *m;
}
#endif
#endif
