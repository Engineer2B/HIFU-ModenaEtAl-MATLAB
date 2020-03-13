pointPlt = @(x,y,z) line([x x],[y y],[z z], 'Marker', '+');

figure;
hold on;
myRay = ray(ri);
line([myRay.start(1) myRay.end(1)],...
	[myRay.start(2) myRay.end(2)],...
	[myRay.start(3) myRay.end(3)]);
for indexLambda = 1:length(lambda_interesting)
	myLambda = lambda_interesting(indexLambda);
	pointPlt(myRay.start(1)+myLambda*myRay.Vray(1),...
		myRay.start(2)+myLambda*myRay.Vray(2),...
		myRay.start(3)+myLambda*myRay.Vray(3));
end
for indexX = 1:length(xxb)
	X = xxb(indexX);
	for indexY = 1:length(yyb)
		Y = yyb(indexY);
		for indexZ = 1:length(zzb)
			Z = zzb(indexZ);
			pointPlt(X,Y,Z);
		end
	end
end
hold off;