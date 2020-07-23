function [BER,SER] = Prevoting_MMSE(n_symbols , n_iterations , bits_all , n_all , H_all , par)

n_errors=zeros(length(par.SNRdB_list),1); %Numero de Erros para cada SNR -> Simbolos
n_errors_bits=zeros(length(par.SNRdB_list),1); %Numero de Erros para cada SNR -> Bits

for loop1 = 1:n_iterations   
    H = H_all(:,:,loop1); %Matriz Canal com Distribuição Normal com media zero e varianza unitaria
    
    for ind_db = 1:length(par.SNRdB_list)
        sigma = sqrt(par.MT*par.Es*10^(-par.SNRdB_list(ind_db)/10));        
        %Inicialização Erro
        err = 0;
        err_bits = 0;
        
        %parfor loop2=1:n_symbols
        for loop2 = 1:n_symbols
            % generate transmit symbol
            idx = bi2de(bits_all(:,:,loop2),'left-msb')+1;
            x = par.symbols(idx).'; 
            % Sinal Recebido = Canal*SinalEnviado + Ruido
            y = H*x + sigma*n_all(:,loop2,loop1);    
            %% Começo Detecção Pre-Voting
            %Pre-Voting
            R = par.MT - par.MR;  %Saber como vou agrupar 
            xp = 0:2^(R*par.Q)-1;
            for ii = 1:size(xp,2) 
                binario = de2bi(xp(:,ii) ,'left-msb',R*par.Q);
                binario = reshape(binario,R,[]);     
                idx = bi2de(binario,'left-msb')+1;
                xp_symbol(:,ii) = par.symbols(idx).';
            end
            Hp = H(:,1:R);               %Canal Prevoting
            Hpxp = Hp*xp_symbol;
            %Canal Posvoting
            Hq = H(:,R+1:end);                
            Gq = pinv( Hq'*Hq + (sigma^2/par.Es)*eye(par.MR) )*Hq'; %Calculo Matriz Inversa de Canal de Postvoting Uso MMSE
            %Começa Algoritmo
            for i = 1:size(Hpxp,2)
                r = y - Hpxp(:,i);
                xq_hat = Gq*r;  %Faço o ZF
                [none,idxhat] = min(abs(xq_hat*ones(1,length(par.symbols))-ones(size(xq_hat,1),1)*par.symbols).^2,[],2);
                xq_hat = par.symbols(idxhat).'; %Decodificação do sq_hat - Quantização                  
                x_hat(:,i) = [xp_symbol(:,i);xq_hat]; %Gero Todos os Candidatos    
            end
            yy = repmat(y,1,size(x_hat,2));
            H_linha = [Hp Hq];
            [non,posicion] =  min(sum(abs(yy - H_linha*x_hat).^2)); %Dist. Euclideana
            x_decod = x_hat(:,posicion); 
            [non,idxhat] = min(abs(x_decod*ones(1,length(par.symbols))-ones(par.MT,1)*par.symbols).^2,[],2);  
            bits_hat = par.bits(idxhat,:);            
            %Conteo Erros
            err = err + sum(x~=x_decod); %Erro De Simbolos
            err_bits = err_bits + sum(sum(bits_all(:,:,loop2)~=bits_hat)); %Erro De Bits
        end %fim Loop2-> n_symbols
        n_errors(ind_db)=err;
        n_errors_bits(ind_db)=err_bits;
    end
end
            
SER = n_errors/(par.MT*n_iterations*n_symbols); %SER         
BER = n_errors_bits/(par.MT*par.Q*n_iterations*n_symbols); %BER         