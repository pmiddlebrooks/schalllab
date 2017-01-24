function env = get_environment


switch matlabroot
    case '/Applications/MATLAB_R2016a.app'
        hostname = char( getHostName( java.net.InetAddress.getLocalHost ) );
        switch hostname
            case 'pmiddleb-Mac-Pro.local'
                env = 'work';
            otherwise
                env = 'home';
        end
    otherwise
        env = 'accre';
        disp('I think you are running on ACCRE- if not, check get_environment.m. Maybe you updated matlab version')
end

