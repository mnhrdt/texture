function x = asc_read(filename)
	f = fopen(filename, "r");
	v = fscanf(f, "%g");
	fclose(f);
	x = reshape(v(5:end), v(1), v(2), v(4));
end
