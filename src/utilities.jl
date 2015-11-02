# utilities functions
using PyPlot

function plot_tps(res_arr,dim)
    tps_val = zeros(8,3)
    for j = 1:3
        for i = 1:8
            tps_val[i,j] = res_arr[j+(i-1)*3]["tps"][dim]
        end
    end
    Nval=collect(250:500:3750)
    etaval=[1.0;1.001;1.01]
    for j = 1:3
        PyPlot.plot(Nval,vec(tps_val[:,j]),"-*")
    end
    legend(etaval,loc="best")
    xlabel("N")
    ylabel("TPS")
end

function plot_ene(res_arr)
    homo = zeros(8,3)
    lumo = zeros(8,3)
    for j = 1:3
        for i = 1:8
            homo[i,j] = res_arr[j+(i-1)*3]["energies"][end/2]
            lumo[i,j] = res_arr[j+(i-1)*3]["energies"][end/2+1]
        end
    end
    Nval=collect(250:500:3750)
    etaval=[1.0;1.001;1.01]
    colors = zeros(3,8)
    for i = 1:8
        colors[:,i] = rand(3)
        plot(etaval, vec(homo[i,:]), color=colors[:,i],"-*")
    end
    legend(Nval,loc="best")
    for i = 1:8
        plot(etaval, vec(lumo[i,:]), color=colors[:,i],"-*")
    end
    xlabel("eta")
    ylabel("E")
end
