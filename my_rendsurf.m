%% start building the figure
h = figure;
set(h, 'color', [1 1 1]);
set(h, 'visible', 'on');
set(h, 'renderer', cfg.renderer);
if ~isempty(cfg.title)
  title(cfg.title);
end

%%% set color and opacity mapping for this figure
if hasfun
  colormap(cfg.funcolormap);
  cfg.funcolormap = colormap;
end
if hasmsk
  cfg.opacitymap = alphamap(cfg.opacitymap);
  alphamap(cfg.opacitymap);
  if ndims(fun)>3 && ndims(msk)==3
    siz = size(fun);
    msk = repmat(msk, [1 1 1 siz(4:end)]);
  end
end

switch cfg.method
  case 'slice'
    % set the defaults for method=slice
    cfg.nslices    = ft_getopt(cfg, 'nslices',    20);
    cfg.slicedim   = ft_getopt(cfg, 'slicedim',   3);
    cfg.slicerange = ft_getopt(cfg, 'slicerange', 'auto');



    % white BG => mskana

    % TODO: HERE THE FUNCTION THAT MAKES TO SLICE DIMENSION ALWAYS THE THIRD DIMENSION, AND ALSO KEEP TRANSFORMATION MATRIX UP TO DATE
    % zoiets
    % if hasana; ana = shiftdim(ana,cfg.slicedim-1); end;
    % if hasfun; fun = shiftdim(fun,cfg.slicedim-1); end;
    % if hasmsk; msk = shiftdim(msk,cfg.slicedim-1); end;

    % ADDED BY JM: allow for slicedim different than 3
    switch cfg.slicedim
      case 1
        dim = dim([2 3 1]);
        if hasana, ana = permute(ana,[2 3 1]); end
        if hasfun, fun = permute(fun,[2 3 1]); end
        if hasmsk, msk = permute(msk,[2 3 1]); end
        cfg.slicedim=3;
      case 2
        dim = dim([3 1 2]);
        if hasana, ana = permute(ana,[3 1 2]); end
        if hasfun, fun = permute(fun,[3 1 2]); end
        if hasmsk, msk = permute(msk,[3 1 2]); end
        cfg.slicedim=3;
      otherwise
        % nothing needed
    end

    %%%%% select slices
    if ~ischar(cfg.slicerange)
      ind_fslice = cfg.slicerange(1);
      ind_lslice = cfg.slicerange(2);
    elseif isequal(cfg.slicerange, 'auto')
      if hasfun % default
        if isfield(functional, 'inside')

          insideMask = false(size(fun));
          insideMask(functional.inside) = true;

          ind_fslice = min(find(max(max(insideMask,[],1),[],2)));
          ind_lslice = max(find(max(max(insideMask,[],1),[],2)));
        else
          ind_fslice = min(find(~isnan(max(max(fun,[],1),[],2))));
          ind_lslice = max(find(~isnan(max(max(fun,[],1),[],2))));
        end
      elseif hasana % if only ana, no fun
        ind_fslice = min(find(max(max(ana,[],1),[],2)));
        ind_lslice = max(find(max(max(ana,[],1),[],2)));
      else
        error('no functional parameter and no anatomical parameter, can not plot');
      end
    else
      error('do not understand cfg.slicerange');
    end
    ind_allslice = linspace(ind_fslice,ind_lslice,cfg.nslices);
    ind_allslice = round(ind_allslice);
    % make new ana, fun, msk, mskana with only the slices that will be plotted (slice dim is always third dimension)
    if hasana; new_ana = ana(:,:,ind_allslice); clear ana; ana=new_ana; clear new_ana; end;
    if hasfun; new_fun = fun(:,:,ind_allslice); clear fun; fun=new_fun; clear new_fun; end;
    if hasmsk; new_msk = msk(:,:,ind_allslice); clear msk; msk=new_msk; clear new_msk; end;
    % if hasmskana; new_mskana = mskana(:,:,ind_allslice); clear mskana; mskana=new_mskana; clear new_mskana; end;

    % update the dimensions of the volume
    if hasana; dim=size(ana); else dim=size(fun); end;

    %%%%% make a "quilt", that contain all slices on 2D patched sheet
    % Number of patches along sides of Quilt (M and N)
    % Size (in voxels) of side of patches of Quilt (m and n)

    % take care of a potential singleton 3rd dimension
    if numel(dim)<3
      dim(end+1:3) = 1;
    end

    %if cfg.slicedim~=3
    %  error('only supported for slicedim=3');
    %end


    m = dim(1);
    n = dim(2);
    M = ceil(sqrt(dim(3)));
    N = ceil(sqrt(dim(3)));
    num_patch = N*M;

    num_slice = (dim(cfg.slicedim));
    num_empt = num_patch-num_slice;
    % put empty slides on ana, fun, msk, mskana to fill Quilt up
    if hasana; ana(:,:,end+1:num_patch)=0; end;
    if hasfun; fun(:,:,end+1:num_patch)=0; end;
    if hasmsk; msk(:,:,end+1:num_patch)=0; end;
    % if hasmskana; mskana(:,:,end:num_patch)=0; end;
    % put the slices in the quilt
    for iSlice = 1:num_slice
      xbeg = floor((iSlice-1)./M);
      ybeg = mod(iSlice-1, M);
      if hasana
        quilt_ana(ybeg*m+1:(ybeg+1)*m, xbeg*n+1:(xbeg+1)*n)=ana(:,:,iSlice);
      end
      if hasfun
        quilt_fun(ybeg*m+1:(ybeg+1)*m, xbeg*n+1:(xbeg+1)*n)=fun(:,:,iSlice);
      end
      if hasmsk
        quilt_msk(ybeg.*m+1:(ybeg+1)*m, xbeg*n+1:(xbeg+1)*n)=msk(:,:,iSlice);
      end
      %     if hasmskana
      %       quilt_mskana(ybeg.*m+1:(ybeg+1).*m, xbeg.*n+1:(xbeg+1).*n)=mskana(:,:,iSlice);
      %     end
    end
    % make vols and scales, containes volumes to be plotted (fun, ana, msk), added by ingnie
    if hasana; vols2D{1} = quilt_ana; scales{1} = []; end; % needed when only plotting ana
    if hasfun; vols2D{2} = quilt_fun; scales{2} = [fcolmin fcolmax]; end;
    if hasmsk; vols2D{3} = quilt_msk; scales{3} = [opacmin opacmax]; end;

    % the transpose is needed for displaying the matrix using the MATLAB image() function
    if hasana;             ana = vols2D{1}'; end;
    if hasfun && ~doimage; fun = vols2D{2}'; end;
    if hasfun &&  doimage; fun = permute(vols2D{2},[2 1 3]); end;
    if hasmsk;             msk = vols2D{3}'; end;

    if hasana
      % scale anatomy between 0 and 1
      fprintf('scaling anatomy\n');
      amin = min(ana(:));
      amax = max(ana(:));
      ana = (ana-amin)./(amax-amin);
      clear amin amax;
      % convert anatomy into RGB values
      ana = cat(3, ana, ana, ana);
      ha = imagesc(ana);
    end
    hold on

    if hasfun

      if doimage
        hf = image(fun);
      else
        hf = imagesc(fun);
        try
          caxis(scales{2});
        end
        % apply the opacity mask to the functional data
        if hasmsk
          % set the opacity
          set(hf, 'AlphaData', msk)
          set(hf, 'AlphaDataMapping', 'scaled')
          try
            alim(scales{3});
          end
        elseif hasana
          set(hf, 'AlphaData', 0.5)
        end

      end
    end

    axis equal
    axis tight
    axis xy
    axis off

    if istrue(cfg.colorbar)
      if hasfun
        % use a normal MATLAB coorbar
        hc = colorbar;
        set(hc, 'YLim', [fcolmin fcolmax]);
      else
        warning('no colorbar possible without functional data')
      end
    end

  case 'ortho'
    % set the defaults for method=ortho
    cfg.location            = ft_getopt(cfg, 'location',            'auto');
    cfg.locationcoordinates = ft_getopt(cfg, 'locationcoordinates', 'head');
    cfg.crosshair           = ft_getopt(cfg, 'crosshair',           'yes');
    cfg.axis                = ft_getopt(cfg, 'axis',                'on');
    cfg.queryrange          = ft_getopt(cfg, 'queryrange',          3);

    if ~ischar(cfg.location)
      if strcmp(cfg.locationcoordinates, 'head')
        % convert the headcoordinates location into voxel coordinates
        loc = inv(functional.transform) * [cfg.location(:); 1];
        loc = round(loc(1:3));
      elseif strcmp(cfg.locationcoordinates, 'voxel')
        % the location is already in voxel coordinates
        loc = round(cfg.location(1:3));
      else
        error('you should specify cfg.locationcoordinates');
      end
    else
      if isequal(cfg.location, 'auto')
        if hasfun
          if isequal(cfg.funcolorlim, 'maxabs');
            loc = 'max';
          elseif isequal(cfg.funcolorlim, 'zeromax');
            loc = 'max';
          elseif isequal(cfg.funcolorlim, 'minzero');
            loc = 'min';
          else % if numerical
            loc = 'max';
          end
        else
          loc = 'center';
        end;
      else
        loc = cfg.location;
      end
    end

    % determine the initial intersection of the cursor (xi yi zi)
    if ischar(loc) && strcmp(loc, 'min')
      if isempty(cfg.funparameter)
        error('cfg.location is min, but no functional parameter specified');
      end
      [dummy, minindx] = min(fun(:));
      [xi, yi, zi] = ind2sub(dim, minindx);
    elseif ischar(loc) && strcmp(loc, 'max')
      if isempty(cfg.funparameter)
        error('cfg.location is max, but no functional parameter specified');
      end
      [dummy, maxindx] = max(fun(:));
      [xi, yi, zi] = ind2sub(dim, maxindx);
    elseif ischar(loc) && strcmp(loc, 'center')
      xi = round(dim(1)/2);
      yi = round(dim(2)/2);
      zi = round(dim(3)/2);
    elseif ~ischar(loc)
      % using nearest instead of round ensures that the position remains within the volume
      xi = nearest(1:dim(1), loc(1));
      yi = nearest(1:dim(2), loc(2));
      zi = nearest(1:dim(3), loc(3));
    end

    xi = round(xi); xi = max(xi, 1); xi = min(xi, dim(1));
    yi = round(yi); yi = max(yi, 1); yi = min(yi, dim(2));
    zi = round(zi); zi = max(zi, 1); zi = min(zi, dim(3));

    % axes settings
    if strcmp(cfg.axisratio, 'voxel')
      % determine the number of voxels to be plotted along each axis
      axlen1 = dim(1);
      axlen2 = dim(2);
      axlen3 = dim(3);
    elseif strcmp(cfg.axisratio, 'data')
      % determine the length of the edges along each axis
      [cp_voxel, cp_head] = cornerpoints(dim, functional.transform);
      axlen1 = norm(cp_head(2,:)-cp_head(1,:));
      axlen2 = norm(cp_head(4,:)-cp_head(1,:));
      axlen3 = norm(cp_head(5,:)-cp_head(1,:));
    elseif strcmp(cfg.axisratio, 'square')
      % the length of the axes should be equal
      axlen1 = 1;
      axlen2 = 1;
      axlen3 = 1;
    end

    % this is the size reserved for subplot h1, h2 and h3
    h1size(1) = 0.82*axlen1/(axlen1 + axlen2);
    h1size(2) = 0.82*axlen3/(axlen2 + axlen3);
    h2size(1) = 0.82*axlen2/(axlen1 + axlen2);
    h2size(2) = 0.82*axlen3/(axlen2 + axlen3);
    h3size(1) = 0.82*axlen1/(axlen1 + axlen2);
    h3size(2) = 0.82*axlen2/(axlen2 + axlen3);

    if strcmp(cfg.voxelratio, 'square')
      voxlen1 = 1;
      voxlen2 = 1;
      voxlen3 = 1;
    elseif strcmp(cfg.voxelratio, 'data')
      % the size of the voxel is scaled with the data
      [cp_voxel, cp_head] = cornerpoints(dim, functional.transform);
      voxlen1 = norm(cp_head(2,:)-cp_head(1,:))/norm(cp_voxel(2,:)-cp_voxel(1,:));
      voxlen2 = norm(cp_head(4,:)-cp_head(1,:))/norm(cp_voxel(4,:)-cp_voxel(1,:));
      voxlen3 = norm(cp_head(5,:)-cp_head(1,:))/norm(cp_voxel(5,:)-cp_voxel(1,:));
    end

    %% the figure is interactive, add callbacks
    set(h, 'windowbuttondownfcn', @cb_buttonpress);
    set(h, 'windowbuttonupfcn',   @cb_buttonrelease);
    set(h, 'windowkeypressfcn',   @cb_keyboard);
    set(h, 'CloseRequestFcn',     @cb_cleanup);

    % ensure that this is done in interactive mode
    set(h, 'renderer', cfg.renderer);

    %% create figure handles

    % axis handles will hold the anatomical functional if present, along with labels etc.
    h1 = axes('position',[0.06                0.06+0.06+h3size(2) h1size(1) h1size(2)]);
    h2 = axes('position',[0.06+0.06+h1size(1) 0.06+0.06+h3size(2) h2size(1) h2size(2)]);
    h3 = axes('position',[0.06                0.06                h3size(1) h3size(2)]);

    set(h1, 'Tag', 'ik', 'Visible', cfg.axis, 'XAxisLocation', 'top');
    set(h2, 'Tag', 'jk', 'Visible', cfg.axis, 'YAxisLocation', 'right'); % after rotating in ft_plot_ortho this becomes top
    set(h3, 'Tag', 'ij', 'Visible', cfg.axis);

    set(h1, 'DataAspectRatio',1./[voxlen1 voxlen2 voxlen3]);
    set(h2, 'DataAspectRatio',1./[voxlen1 voxlen2 voxlen3]);
    set(h3, 'DataAspectRatio',1./[voxlen1 voxlen2 voxlen3]);

    % create structure to be passed to gui
    opt               = [];
    opt.dim           = dim;
    opt.ijk           = [xi yi zi];
    opt.h1size        = h1size;
    opt.h2size        = h2size;
    opt.h3size        = h3size;
    opt.handlesaxes   = [h1 h2 h3];
    opt.handlesfigure = h;
    opt.axis          = cfg.axis;
    if hasatlas
      opt.atlas = atlas;
    end
    if hasana
      opt.ana = ana;
    end
    if hasfun
      opt.fun = fun;
    end
    if hasmsk
      opt.msk = msk;
    end
    opt.update        = [1 1 1];
    opt.init          = true;
    opt.usedim        = (isUnstructuredFun==false);
    opt.usepos        = (isUnstructuredFun==true);
    opt.hasatlas      = hasatlas;
    opt.hasfreq       = hasfreq;
    opt.hastime       = hastime;
    opt.hasmsk        = hasmsk;
    opt.hasfun        = hasfun;
    opt.hasana        = hasana;
    opt.qi            = qi;
    opt.tag           = 'ik';
    opt.functional    = functional;
    opt.fcolmin       = fcolmin;
    opt.fcolmax       = fcolmax;
    opt.opacmin       = opacmin;
    opt.opacmax       = opacmax;
    opt.clim          = []; % contrast limits for the anatomy, see ft_volumenormalise
    opt.colorbar      = cfg.colorbar;
    opt.queryrange    = cfg.queryrange;
    opt.funcolormap   = cfg.funcolormap;
    opt.crosshair     = istrue(cfg.crosshair);

    %% do the actual plotting
    setappdata(h, 'opt', opt);
    cb_redraw(h);

    fprintf('\n');
    fprintf('click left mouse button to reposition the cursor\n');
    fprintf('click and hold right mouse button to update the position while moving the mouse\n');
    fprintf('use the arrowkeys to navigate in the current axis\n');


  case 'surface'
    % set the defaults for method=surface
    cfg.downsample     = ft_getopt(cfg, 'downsample',     1);
    cfg.surfdownsample = ft_getopt(cfg, 'surfdownsample', 1);
    cfg.surffile       = ft_getopt(cfg, 'surffile', 'surface_white_both.mat'); % use a triangulation that corresponds with the collin27 anatomical template in MNI coordinates
    cfg.surfinflated   = ft_getopt(cfg, 'surfinflated',  []);
    cfg.sphereradius   = ft_getopt(cfg, 'sphereradius',  []);
    cfg.projvec        = ft_getopt(cfg, 'projvec',       1);
    cfg.projweight     = ft_getopt(cfg, 'projweight',    ones(size(cfg.projvec)));
    cfg.projcomb       = ft_getopt(cfg, 'projcomb',      'mean'); % or max
    cfg.projthresh     = ft_getopt(cfg, 'projthresh',    []);
    cfg.projmethod     = ft_getopt(cfg, 'projmethod',    'nearest');
    cfg.distmat        = ft_getopt(cfg, 'distmat',       []);
    cfg.camlight       = ft_getopt(cfg, 'camlight',      'yes');

    % determine whether the source functional already contains a triangulation
    interpolate2surf = 0;
    if ~isUnstructuredFun
      % no triangulation present: interpolation should be performed
      fprintf('The source functional is defined on a 3D grid, interpolation to a surface mesh will be performed\n');
      interpolate2surf = 1;
    elseif isUnstructuredFun && isfield(functional, 'tri')
      fprintf('The source functional is defined on a triangulated surface, using the surface mesh description in the functional\n');
    elseif isUnstructuredFun
      % add a transform field to the functional
      fprintf('The source functional does not contain a triangulated surface, we may need to interpolate to a surface mesh\n');
      functional.transform = pos2transform(functional.pos);
      interpolate2surf = 1;
    end

    if interpolate2surf,
      % deal with the interpolation
      % FIXME this should be dealt with by ft_sourceinterpolate

      % read the triangulated cortical surface from file
      surf = ft_read_headshape(cfg.surffile);

      if isfield(surf, 'transform'),
        % compute the surface vertices in head coordinates
        surf.pos = ft_warp_apply(surf.transform, surf.pos);
      end

      % downsample the cortical surface
      if cfg.surfdownsample > 1
        if ~isempty(cfg.surfinflated)
          error('downsampling the surface is not possible in combination with an inflated surface');
        end
        fprintf('downsampling surface from %d vertices\n', size(surf.pos,1));
        [temp.tri, temp.pos] = reducepatch(surf.tri, surf.pos, 1/cfg.surfdownsample);
        % find indices of retained patch faces
        [dummy, idx] = ismember(temp.pos, surf.pos, 'rows');
        idx(idx==0)  = [];
        surf.tri = temp.tri;
        surf.pos = temp.pos;
        clear temp
        % downsample other fields
        if isfield(surf, 'curv'),       surf.curv       = surf.curv(idx);       end
        if isfield(surf, 'sulc'),       surf.sulc       = surf.sulc(idx);       end
        if isfield(surf, 'hemisphere'), surf.hemisphere = surf.hemisphere(idx); end
      end

      % these are required
      if ~isfield(functional, 'inside')
        functional.inside = true(dim);
      end

      fprintf('%d voxels in functional data\n', prod(dim));
      fprintf('%d vertices in cortical surface\n', size(surf.pos,1));

      tmpcfg = [];
      tmpcfg.parameter = {cfg.funparameter};
      if ~isempty(cfg.maskparameter)
        tmpcfg.parameter = [tmpcfg.parameter {cfg.maskparameter}];
        maskparameter    = cfg.maskparameter;
      else
        tmpcfg.parameter = [tmpcfg.parameter {'mask'}];
        functional.mask  = msk;
        maskparameter    = 'mask'; % temporary variable
      end
      tmpcfg.interpmethod = cfg.projmethod;
      tmpcfg.distmat      = cfg.distmat;
      tmpcfg.sphereradius = cfg.sphereradius;
      tmpcfg.projvec      = cfg.projvec;
      tmpcfg.projcomb     = cfg.projcomb;
      tmpcfg.projweight   = cfg.projweight;
      tmpcfg.projthresh   = cfg.projthresh;
      tmpdata             = ft_sourceinterpolate(tmpcfg, functional, surf);

      if hasfun, val      = getsubfield(tmpdata, cfg.funparameter);  val     = val(:);     end
      if hasmsk, maskval  = getsubfield(tmpdata, maskparameter);     maskval = maskval(:); end

      if ~isempty(cfg.projthresh),
        maskval(abs(val) < cfg.projthresh*max(abs(val(:)))) = 0;
      end

    else
      surf     = [];
      surf.pos = functional.pos;
      surf.tri = functional.tri;

      % if hasfun, val     = fun(functional.inside(:)); end
      % if hasmsk, maskval = msk(functional.inside(:)); end
      if hasfun, val     = fun(:); end
      if hasmsk, maskval = msk(:); end

    end

    if ~isempty(cfg.surfinflated)
      if ~isstruct(cfg.surfinflated)
        % read the inflated triangulated cortical surface from file
        surf = ft_read_headshape(cfg.surfinflated);
      else
        surf = cfg.surfinflated;
        if isfield(surf, 'transform'),
          % compute the surface vertices in head coordinates
          surf.pos = ft_warp_apply(surf.transform, surf.pos);
        end
      end
    end

    %------do the plotting
    cortex_light = [0.781 0.762 0.664];
    cortex_dark  = [0.781 0.762 0.664]/2;
    if isfield(surf, 'curv')
      % the curvature determines the color of gyri and sulci
      color = surf.curv(:) * cortex_dark + (1-surf.curv(:)) * cortex_light;
    else
      color = repmat(cortex_light, size(surf.pos,1), 1);
    end

    h1 = patch('Vertices', surf.pos, 'Faces', surf.tri, 'FaceVertexCData', color, 'FaceColor', 'interp');
    set(h1, 'EdgeColor', 'none');
    axis   off;
    axis vis3d;
    axis equal;

    if hasfun
      h2 = patch('Vertices', surf.pos, 'Faces', surf.tri, 'FaceVertexCData', val, 'FaceColor', 'interp');
      set(h2, 'EdgeColor', 'none');
      try
        caxis(gca,[fcolmin fcolmax]);
      end
      colormap(cfg.funcolormap);
      if hasmsk
        set(h2, 'FaceVertexAlphaData', maskval);
        set(h2, 'FaceAlpha',          'interp');
        set(h2, 'AlphaDataMapping',   'scaled');
        try
          alim(gca, [opacmin opacmax]);
        end
        alphamap(cfg.opacitymap);
      end
    end

    lighting gouraud

    if istrue(cfg.camlight)
      camlight
    end

    if istrue(cfg.colorbar)
      if hasfun
        % use a normal MATLAB colorbar
        hc = colorbar;
        set(hc, 'YLim', [fcolmin fcolmax]);
      else
        warning('no colorbar possible without functional data')
      end
    end
