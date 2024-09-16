function graphics_setup(plot_type)

switch plot_type
case '1by1'
set(gcf, 'units', 'centimeters', 'Position',  [12, 3, 7, 5.5]);
set(gca,'fontsize',7)

case '1by1.5'
set(gcf, 'units', 'centimeters', 'Position',  [12, 3, 8.5, 5.5]);
set(gca,'fontsize',7)

case '1by2'
set(gcf, 'units', 'centimeters', 'Position',  [3, 3, 11.5, 5]);
set(gca,'fontsize',7)

case '1by2b'
set(gcf, 'units', 'centimeters', 'Position',  [3, 3, 11, 5]);
set(gca,'fontsize',7)

case '2by3'
set(gcf, 'units', 'centimeters', 'Position',  [3, 3, 17, 11]);
set(gca,'fontsize',7)

case '2by2'
set(gcf, 'units', 'centimeters', 'Position',  [12, 3, 11.5, 10.5]);
set(gca,'fontsize',7)

case '3by4'
set(gcf, 'units', 'centimeters', 'Position',  [3, 3, 18, 15]);
set(gca,'fontsize',5)

end

return