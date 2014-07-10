#include "time_expr.hpp"

#include <limits>
#include <deque>
#include <boost/algorithm/string/trim.hpp>

using boost::lexical_cast;
using boost::trim_copy;

namespace realevol {

// Singleton pool of mu::Parser objects
class muParser_pool {
    std::deque<mu::Parser*> parsers;
    static constexpr const std::size_t chunk_size = 5;

public:
    mu::Parser* acquire()
    {
        if(!parsers.size()){
            for(std::size_t n=0; n<chunk_size; ++n) parsers.push_back(new mu::Parser);
        }
        auto p = parsers.back();
        parsers.pop_back();
        return p;
    }

    void release(mu::Parser* p)
    {
        if(p!=nullptr) parsers.push_back(p);
    }
} pool;

time_expr::time_expr(std::string const& expr) :
    parser(pool.acquire())
{
    parser->DefineVar("t", &arg);
    parser->SetExpr(trim_copy(expr));
    if(!parser->GetUsedVar().size()){
        arg = parser->Eval();
        pool.release(parser);
        parser = nullptr;
    }
}

time_expr::time_expr(const char* expr) : time_expr(std::string(expr))
{}

time_expr::time_expr(time_expr::result_type r) :
    arg(r), parser(nullptr)
{}

time_expr::time_expr() : time_expr(result_type(0))
{}

time_expr::time_expr(time_expr const& te) :
    arg(te.arg), parser(nullptr)
{
    if(te.parser != nullptr){
        parser = pool.acquire();
        parser->DefineVar("t",&arg);
        parser->SetExpr(trim_copy(te.parser->GetExpr()));
    }
}

time_expr::~time_expr()
{
    pool.release(parser);
}

auto time_expr::operator()(double t) const -> result_type
{
    if(is_constant(*this))
        return arg;
    else {
        arg = t;
        return parser->Eval();
    }
}

std::ostream& operator<<(std::ostream& os, time_expr const& te)
{
    if(is_constant(te))
        os << te.arg;
    else
        os << trim_copy(te.parser->GetExpr());
    
    return os;
}

auto time_expr::operator-() const -> time_expr
{
    if(is_constant(*this))
        return time_expr(-arg);
    else
        return time_expr("-(" + trim_copy(parser->GetExpr()) + ")");
}

auto time_expr::operator=(const time_expr& te) -> time_expr
{
    if(is_constant(te)){
        arg = te.arg;
        pool.release(parser);
        parser = nullptr;
    }else{
        if(is_constant(*this)){
            parser = pool.acquire();
            parser->DefineVar("t",&arg);            
        }
        parser->SetExpr(trim_copy(te.parser->GetExpr()));
    }
    return *this;
}

auto time_expr::operator=(std::string const& expr) -> time_expr
{
    if(is_constant(*this))
    {
        parser = pool.acquire();
    }
    parser->DefineVar("t", &arg);
    parser->SetExpr(trim_copy(expr));
    if(!parser->GetUsedVar().size()){
        arg = parser->Eval();
        pool.release(parser);
        parser = nullptr;
    }
    
    return *this;
}

auto time_expr::operator=(const char* expr) -> time_expr
{
    *this = std::string(expr);
    
    return *this;
}

auto time_expr::operator=(result_type r) -> time_expr
{
    if(!is_constant(*this)){
        pool.release(parser);
        parser = nullptr;
    }
    arg = r;
    
    return *this;
}

auto time_expr::operator+=(const time_expr& te) -> time_expr
{
    bool is_const = is_constant(*this);
    if(is_const && is_constant(te))
        arg += te.arg;
    else {
        if(is_const){
            parser = pool.acquire();
            parser->DefineVar("t",&arg);
        }
        parser->SetExpr(
            "(" +
            (is_const ? lexical_cast<std::string>(arg) : trim_copy(parser->GetExpr())) +
            ")+(" +
            (is_constant(te) ? lexical_cast<std::string>(te.arg) : trim_copy(te.parser->GetExpr())) +
            ")"
        );
    }
    
    return *this;
}

auto time_expr::operator-=(const time_expr& te) -> time_expr
{
    bool is_const = is_constant(*this);
    if(is_const && is_constant(te))
        arg -= te.arg;
    else {
        if(is_const){
            parser = pool.acquire();
            parser->DefineVar("t",&arg);
        }
        parser->SetExpr(
            "(" +
            (is_const ? lexical_cast<std::string>(arg) : trim_copy(parser->GetExpr())) +
            ")-(" +
            (is_constant(te) ? lexical_cast<std::string>(te.arg) : trim_copy(te.parser->GetExpr())) +
            ")"
        );
    }
    
    return *this;
}

auto time_expr::operator*=(const time_expr& te) -> time_expr
{
    bool is_const = is_constant(*this);
    if(is_const && is_constant(te))
        arg *= te.arg;
    else {
        if(is_const){
            parser = pool.acquire();
            parser->DefineVar("t",&arg);
        }
        parser->SetExpr(
            "(" +
            (is_const ? lexical_cast<std::string>(arg) : trim_copy(parser->GetExpr())) +
            ")*(" +
            (is_constant(te) ? lexical_cast<std::string>(te.arg) : trim_copy(te.parser->GetExpr())) +
            ")"
        );
    }
    
    return *this;
}

auto time_expr::operator/=(const time_expr& te) -> time_expr
{
    bool is_const = is_constant(*this);
    if(is_const && is_constant(te))
        arg /= te.arg;
    else {
        if(is_const){
            parser = pool.acquire();
            parser->DefineVar("t",&arg);
        }
        parser->SetExpr(
            "(" +
            (is_const ? lexical_cast<std::string>(arg) : trim_copy(parser->GetExpr())) +
            ")/(" +
            (is_constant(te) ? lexical_cast<std::string>(te.arg) : trim_copy(te.parser->GetExpr())) +
            ")"
        );
    }
    
    return *this;
}

bool time_expr::operator==(const time_expr& te) const
{
    bool is_const = is_constant(*this);
    if(is_const != is_constant(te)) return false;
    if(is_const)
        return is_zero(arg-te.arg);
    else
        return parser->GetExpr() == te.parser->GetExpr();
}
    
}