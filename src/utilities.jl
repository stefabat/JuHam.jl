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
        #PyPlot.plot(Nval,vec(tps_val[:,j]),"-*")
        plot(Nval,vec(tps_val[:,j]),"-*")
    end
    legend(etaval,loc="best")
    xlabel("N")
    ylabel("TPS")
end

function plot_polarizability(res_arr,dim)
    pol_val = zeros(8,3)
    for j = 1:3
        for i = 1:8
            pol_val[i,j] = res_arr[j+(i-1)*3]["polarizability"][dim]
        end
    end
    Nval=collect(250:500:3750)
    etaval=[1.0;1.001;1.01]
    for j = 1:3
        #PyPlot.plot(Nval,vec(pol_val[:,j]),"-*")
        plot(Nval,vec(pol_val[:,j]),"-*")
    end
    legend(etaval,loc="best")
    xlabel("N")
    ylabel("Polarizability")
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

# Simple helper function to visualize any topology
function plot_nanotube(topology)
    x     = topology.coords[:,1]
    y     = topology.coords[:,2]
    z     = topology.coords[:,3]
    bonds = topology.bonds
    PyPlot.hold(true)
    for i in bonds
        plot3D([x[i[1]],x[i[2]]],[y[i[1]],y[i[2]]],[z[i[1]],z[i[2]]],"-*k")
    end
    axis(xmin=minimum(x)-1,xmax=maximum(x)+1,ymin=minimum(y)-1,ymax=maximum(y)+1,zmin=minimum(z)-1,zmax=maximum(z)+1,aspect="equal")
    PyPlot.hold(false)
end

function plot_graphene(topology)
    x     = topology.coords[:,1]
    y     = topology.coords[:,2]
    PyPlot.hold(true)
    for i in topology.bonds
        plot([x[i[1]],x[i[2]]],[y[i[1]],y[i[2]]],"-*k")
    end
    axis(xmin=minimum(x)-1,xmax=maximum(x)+1,ymin=minimum(y)-1,ymax=maximum(y)+1,aspect="equal")
    PyPlot.hold(false)
end


