%Test map - decided to stick with plot from outputMetadata... function -
%couldn't get this to work any better.


p = projcrs(4087);
load coastlines coastlon coastlat
land = readgeotable("landareas.shp");

shape1 = geopolyshape([90, 90, 50, 50, 40, 40, 50, 50, 90], [-180, 180, 180, 10, 10, -80, -80, -180, -180]);
shape2 = geopolyshape([-50, -50, -90, -90, -50], [-180 180 180 -180 -180])

figure;


%Set up plot of continents
newmap(p);
geoplot(land, "FaceColor",[0.5 0.5 0.5],"EdgeColor",'k')
hold on

%Plot regions that are being excluded
geoplot(shape1, "FaceColor", "k", "FaceAlpha", 0.5)
geoplot(shape2, "FaceColor", "k", "FaceAlpha", 0.5)

geoscatter(lat1, lon1, 'MarkerEdgeColor', subset1Color, 'LineWidth', 1, 'Marker','square')
geoscatter(lat2, lon2, 'MarkerEdgeColor', subset2Color, 'LineWidth', 1, 'Marker','square')
