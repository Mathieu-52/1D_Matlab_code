function [] = my_plotter(x_c,phi,c1,c2,c_one_scalar,c_one_scalar_i,srf,srf_os,T,T_i,T_1,T_2,my_plot)
plot(x_c,1*(phi),'-k','linewidth',my_plot.line_width);
hold on
if(my_plot.c1)
    plot(x_c,c1,'-.b','linewidth',my_plot.line_width);
    hold on
end
if(my_plot.c2)
    plot(x_c,c2,':r','linewidth',my_plot.line_width);
    hold on
end
if(my_plot.c_one_scalar)
    plot(x_c,c_one_scalar,'-.r','linewidth',my_plot.line_width,'markersize',my_plot.marker_size);
    hold on
end
if(my_plot.c_one_scalar_i)
    plot(x_c,c_one_scalar_i,':m','linewidth',my_plot.line_width,'markersize',my_plot.marker_size);
    hold on
end
if(my_plot.c1pc2)
    plot(x_c,c1+c2,'--m','linewidth',my_plot.line_width,'markersize',my_plot.marker_size);
    hold on
end
if(my_plot.srf)
    plot(x_c,srf,'-m','linewidth',my_plot.line_width,'markersize',my_plot.marker_size);
    hold on
end
if(my_plot.srf_total)
    plot(x_c,srf+c1+c2,':k','linewidth',my_plot.line_width,'markersize',my_plot.marker_size);
    hold on
end
if(my_plot.srf_os)
    plot(x_c,srf_os,':g','linewidth',my_plot.line_width,'markersize',my_plot.marker_size);
    hold on
end
if(my_plot.T)
    p_l=plot(x_c,T,'--k','linewidth',my_plot.line_width,'markersize',my_plot.marker_size);
    hold on
end
if(my_plot.T_i)
    p_l=plot(x_c,T_i,'--k','linewidth',my_plot.line_width,'markersize',my_plot.marker_size);
    hold on
end
if(my_plot.T_1)
    p_l=plot(x_c,T_1,'--k','linewidth',my_plot.line_width,'markersize',my_plot.marker_size);
    hold on
end
if(my_plot.T_2)
    p_l=plot(x_c,T_2,':r','linewidth',my_plot.line_width,'markersize',my_plot.marker_size);
    hold on
end
end