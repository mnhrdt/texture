function asc_write(filename, x)
	f = fopen(filename, "w");
	[w, h, pd] = size(x);
	fprintf(f, "%g %g 1 %g\n", w, h, pd);
	fprintf(f, "%g\n", x);
	fclose(f);
end
