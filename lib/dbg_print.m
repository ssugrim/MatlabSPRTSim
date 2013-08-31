function dbg_print(x, DEBUG)
%DBG_PRINT prints only if DEBUG = true
%   Similar to logger.debug, x is assumed to be an str
if DEBUG == true
    sprintf(x)
end
    

end

