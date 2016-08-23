classdef oriented_prism
    %ORIENTED_PRISM Prism whose orientation is not the standard basis
    
    properties
        corner, diagonal, xh, yh, zh
    end
    
    methods
        function obj = oriented_prism(corner, diagonal, xh, yh, zh)
            obj.corner = corner;
            obj.diagonal = diagonal;
            obj.xh = xh;
            obj.yh = yh;
            obj.zh = zh;
        end
        
        function gz = oriented_prism_gz(self, eval_pt)
            %ORIENTED_PRISM_GZ GZ calculation for an arbitrarily oriented right
            %rectangular prizm
            %   eval_pt in standard basis
            
            % Calculate the full gravitational acceleration in the euclidean frame of
            % the prism
            
            % Find the eval_pt location in the prism's euclidean orientation
            offset_pt = eval_pt - self.corner;
            
            A = [self.xh, self.yh];
            xy = pinv(A) * offset_pt; % Moore Penrose pseudoinverse
            z = dot(offset_pt - A * xy, self.zh);
            
            corner_offset = -[xy; z];
            
            bounds = [corner_offset, corner_offset + self.diagonal];
            
            g = full_g_vector(bounds(1, :), bounds(2, :), bounds(3, :));
            z_dir = [0; 0; 1];
            z_dir_in_prism_basis = [dot(z_dir, self.xh); dot(z_dir, self.yh); dot(z_dir, self.zh)];
            gz = dot(g, z_dir_in_prism_basis);
        end
        
        function render(self)
            render_prism(self.corner, self.diagonal, self.xh, self.yh, self.zh);
        end
    end
end

function g = full_g_vector(x, y, z)
g = [gz_mixed(y, z, x);
     gz_mixed(z, x, y);
     gz_mixed(x, y, z)];
end

function gz_val = gz_mixed(xb, yb, zb)
    gz_val = gz(xb(1), xb(2), yb(1), yb(2), zb(1), zb(2));
end
