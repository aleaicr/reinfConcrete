function [] = createCool3DFigure(ID_data)
% This function create a plot3D figure using M, curvature and c

M = ID_data.M;
curvature = ID_data.curvature;
c = ID_data.c;
M_neg = ID_data.M_neg;
curvature_neg = ID_data.curvature_neg;
c_neg = ID_data.c_neg ;
N_partitions = size(M,2);
Section = ID_data.Section;
N = Section.Pu_;
N_vect = linspace(min(N), max(N), N_partitions)*1000; % kgf

%% Plot
colorPalette = colormap(winter(N_partitions));
colorPalette = colorPalette / max(colorPalette(:));
figure1 = figure('Color',[1 1 1],'position',[680 341 968 637]);
axes1 = axes('Parent',figure1);
hold on
legends = cell(2*N_partitions, 1);
for i = 1:N_partitions
    color_alea = colorPalette(i,:);
    plot3([0; curvature(:,i)], [0; c(:,i)], [0; M(:,i)], 'linewidth', 3, 'Color', color_alea)
    plot3([0; curvature_neg(:,i)], [0; c_neg(:,i)], [0; M_neg(:,i)],'linewidth', 3, 'Color', color_alea)
    legends{2*i-1} = ['N = ' num2str(N_vect(i)/1000) ' [tonf]'];
    legends{2*i} = '';
end
hold off
xlabel('Curvature (phi) [1/cm]')
zlabel('Moment (M) [tonf]')
ylabel('Neutral Axis Depth (c) [cm]')
grid on
set(axes1, 'FontSize', 20);
legend(cellstr(legends));

end

