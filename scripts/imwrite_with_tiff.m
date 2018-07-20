function imwrite_with_tiff(img, filename)
     iio_write(filename, img);
%[~, ~, ext] = fileparts(filename);
%if (strcmp(ext, '.tiff') || strcmp(ext, '.tif'))
%    t = Tiff(filename, 'w');
%    tagstruct.ImageLength = size(img, 1);
%    tagstruct.ImageWidth = size(img, 2);
%    tagstruct.Compression = Tiff.Compression.None;
%    %tagstruct.Compression = Tiff.Compression.LZW;        % compressed
%    tagstruct.SampleFormat = Tiff.SampleFormat.IEEEFP;
%    tagstruct.Photometric = Tiff.Photometric.MinIsBlack;
%    tagstruct.BitsPerSample =  64;                        % double data
%    tagstruct.SamplesPerPixel = size(img,3);
%    tagstruct.PlanarConfiguration = Tiff.PlanarConfiguration.Chunky;
%    t.setTag(tagstruct);
%    t.write((img));
%    t.close();
%else
%    imwrite(img, filename);
end

