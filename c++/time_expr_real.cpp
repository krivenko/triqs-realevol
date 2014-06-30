#include "time_expr.hpp"

#include <limits>
#include <deque>
#include <boost/algorithm/string/trim.hpp>

using boost::lexical_cast;
using boost::trim_copy;

namespace realevol {

// Singleton pool of mu::Parser objects
class muParser_pool {
    typedef std::deque<mu::Parser*> pool_t;
    static pool_t pool;
    static const std::size_t chunk_size = 5;

public:
    static mu::Parser* acquire()
    {
        if(!pool.size()){
            for(std::size_t n=0; n<chunk_size; ++n) pool.push_back(new mu::Parser);
        }
        auto p = pool.back();
        pool.pop_back();
        return p;
    }

    static void release(mu::Parser* p)
    {
        if(p!=nullptr) pool.push_back(p);
    }
};
muParser_pool::pool_t muParser_pool::pool = muParser_pool::pool_t();

time_expr<false>::time_expr(std::string const& expr) :
    parser(muParser_pool::acquire())
{
    parser->DefineVar("t", &arg);
    parser->SetExpr(trim_copy(expr));
    if(!parser->GetUsedVar().size()){
        arg = parser->Eval();
        muParser_pool::release(parser);
        parser = nullptr;
    }
}

time_expr<false>::time_expr(const char* expr) :
    parser(muParser_pool::acquire())
{
    parser->DefineVar("t", &arg);
    parser->SetExpr(expr);
    if(!parser->GetUsedVar().size()){
        arg = parser->Eval();
        muParser_pool::release(parser);
        parser = nullptr;
     }
}

time_expr<false>::time_expr(time_expr::result_type r) :
    arg(r), parser(nullptr)
{}

time_expr<false>::time_expr(int i) :
    arg(i), parser(nullptr)
{}

time_expr<false>::time_expr() :
    arg(0), parser(nullptr)
{}

time_expr<false>::time_expr(time_expr const& te) :
    arg(te.arg), parser(nullptr)
{
    if(te.parser != nullptr){
        parser = muParser_pool::acquire();
        parser->DefineVar("t",&arg);
        parser->SetExpr(trim_copy(te.parser->GetExpr()));
    }
}

time_expr<false>::~time_expr()
{
    muParser_pool::release(parser);
}

bool time_expr<false>::is_constant() const
{
    return parser==nullptr;
}

bool time_expr<false>::is_zero() const
{
    return is_constant() && std::fabs(arg) < comparison_tolerance;
}

auto time_expr<false>::operator()(double t) const -> result_type
{
    if(is_constant())
        return arg;
    else {
        arg = t;
        return parser->Eval();
    }
}

std::ostream& operator<<(std::ostream& os, time_expr<false> const& te)
{
    if(te.is_constant())
        os << te.arg;
    else
        os << trim_copy(te.parser->GetExpr());
    
    return os;
}

auto time_expr<false>::operator-() -> time_expr
{
    if(is_constant())
        return time_expr(-arg);
    else
        return time_expr("-(" + parser->GetExpr() + ")");
}

auto time_expr<false>::operator=(const time_expr& te) -> time_expr
{
    if(te.is_constant()){
        arg = te.arg;
        muParser_pool::release(parser);
        parser = nullptr;
    }else{
        if(is_constant()){
            parser = muParser_pool::acquire();
            parser->DefineVar("t",&arg);            
        }
        parser->SetExpr(trim_copy(te.parser->GetExpr()));
    }
    return *this;
}

auto time_expr<false>::operator+=(const time_expr& te) -> time_expr
{
    bool is_const = is_constant();
    if(is_const && te.is_constant())
        arg += te.arg;
    else {
        if(is_const){
            parser = muParser_pool::acquire();
            parser->DefineVar("t",&arg);
        }
        parser->SetExpr(
            "(" +
            (is_const ? lexical_cast<std::string>(arg) : trim_copy(parser->GetExpr())) +
            ")+(" +
            (te.is_constant() ? lexical_cast<std::string>(te.arg) : trim_copy(te.parser->GetExpr())) +
            ")"
        );
    }
    
    return *this;
}

auto time_expr<false>::operator-=(const time_expr& te) -> time_expr
{
    bool is_const = is_constant();
    if(is_const && te.is_constant())
        arg -= te.arg;
    else {
        if(is_const){
            parser = muParser_pool::acquire();
            parser->DefineVar("t",&arg);
        }
        parser->SetExpr(
            "(" +
            (is_const ? lexical_cast<std::string>(arg) : trim_copy(parser->GetExpr())) +
            ")-(" +
            (te.is_constant() ? lexical_cast<std::string>(te.arg) : trim_copy(te.parser->GetExpr())) +
            ")"
        );
    }
    
    return *this;
}

auto time_expr<false>::operator*=(const time_expr& te) -> time_expr
{
    bool is_const = is_constant();
    if(is_const && te.is_constant())
        arg *= te.arg;
    else {
        if(is_const){
            parser = muParser_pool::acquire();
            parser->DefineVar("t",&arg);
        }
        parser->SetExpr(
            "(" +
            (is_const ? lexical_cast<std::string>(arg) : trim_copy(parser->GetExpr())) +
            ")*(" +
            (te.is_constant() ? lexical_cast<std::string>(te.arg) : trim_copy(te.parser->GetExpr())) +
            ")"
        );
    }
    
    return *this;
}

auto time_expr<false>::operator/=(const time_expr& te) -> time_expr
{
    bool is_const = is_constant();
    if(is_const && te.is_constant())
        arg /= te.arg;
    else {
        if(is_const){
            parser = muParser_pool::acquire();
            parser->DefineVar("t",&arg);
        }
        parser->SetExpr(
            "(" +
            (is_const ? lexical_cast<std::string>(arg) : trim_copy(parser->GetExpr())) +
            ")/(" +
            (te.is_constant() ? lexical_cast<std::string>(te.arg) : trim_copy(te.parser->GetExpr())) +
            ")"
        );
    }
    
    return *this;
}

bool time_expr<false>::operator==(const time_expr& te) const
{
    bool is_const = is_constant();
    if(is_const != te.is_constant()) return false;
    if(is_const)
        return arg == te.arg;
    else
        return parser->GetExpr() == te.parser->GetExpr();
}
    
}