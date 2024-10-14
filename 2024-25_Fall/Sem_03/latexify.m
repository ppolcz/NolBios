function ret = latexify(f,p)
    arguments
        f, p = []
    end

    Eq = "$\left\{\begin{array}{l} \dot x = " + latex(f(1)) + "\\ \dot y = " + latex(f(2)) + " \end{array}\right.$";

    if isempty(p)
        ret = @(~) Eq;
        return
    end

    pnames = cellfun(@(s) {[', $' latex(s) ' = ']}, num2cell(p));

    function ret = lambdaFun(p_val)
        ret = Eq + strjoin(cellfun(@(s,v) string(s) + num2str(v) + "$",pnames,num2cell(p_val)),'');
    end

    ret = @lambdaFun;
end