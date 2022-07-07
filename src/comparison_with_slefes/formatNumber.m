function result=formatNumber(number, digits)
    tmp=num2str(number,2);
    if(contains(tmp,'e'))
        %Number is in engineering format (e.g., 2e+5)
        result=strrep(tmp,'+0',''); %strrep will change 1e+02--> 1e2 (for shorter notation)
    else
        result=num2str(number,['%.',num2str(digits),'f']);
    end
end