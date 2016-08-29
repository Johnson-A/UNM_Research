classdef oriented_prism
    %ORIENTED_PRISM Prism whose orientation is not the standard basis
    %   Corner - standard basis TODO
    
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
        
        function gz = eval_gz_at(self, eval_pt)
            %ORIENTED_PRISM_GZ GZ calculation for an arbitrarily oriented right
            %rectangular prizm
            %   eval_pt in standard basis
            
            % Calculate the full gravitational acceleration in the euclidean frame of
            % the prism
            
            % Find the eval_pt location in the prism's euclidean orientation
            offset_pt = eval_pt - self.corner;
            
            M = [self.xh, self.yh, self.zh];
            corner_offset = -M \ offset_pt;
            z_dir = M \ [0; 0; 1];
            
            bounds = [corner_offset, corner_offset + self.diagonal];
            
            g = full_g_vector(bounds(1, :), bounds(2, :), bounds(3, :));

            gz = dot(g, z_dir);
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
