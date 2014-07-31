function amr_error( status )
%print an error message depending on status
    if (status.value ~= 0)
        status.value
    end
end

