CDF       
      band_lw       band_sw       coeff      	         comment      �This file provides a parameterization of ice particle scattering in the longwave and shortwave RRTM bands,
using the ice particle scattering database of Baran et al. (J. Climate, 2014, 27, 7725–7752),
but with a different functional form for the parameterization. Effective radius is not used.
If ice mass mixing ratio is qi (in kg/kg) and the 1-based 9-element array of coefficients is p, then:
mass extinction coefficient (m2/kg) = qi * {p[1] + p[2] / (1 + qi*p[3])},
single scattering albedo = p[4] + p[5] / (1 + qi*p[6]), and
asymmetry factor = p[7] + p[8] / (1 + qi*p[9]),
where mass extinction coefficient is the total extinction cross section per unit mass of cloudy air.         wavenumber1_lw                  	long_name         (Lower bound wavenumber for longwave band   units         cm-1      @     wavenumber2_lw                  	long_name         (Upper bound wavenumber for longwave band   units         cm-1      @  P   wavenumber1_sw                 	long_name         )Lower bound wavenumber for shortwave band      units         cm-1      8  �   wavenumber2_sw                 	long_name         )Upper bound wavenumber for shortwave band      units         cm-1      8  �   coeff_lw                   	long_name         (Longwave droplet scattering coefficients     @      coeff_sw                  	long_name         )Shortwave droplet scattering coefficients        �  	@A   C�  C�  D� D/  DM  Du  D�  D�� D�� D�  D�  E  E� E� E"� C�  C�  D� D/  DM  Du  D�  D�� D�� D�  D�  E  E� E� E"� EK  E"� EK  Ez  E�P E�� E�0 E� E�� FH� Fz  F�� F� Gp DM  EK  Ez  E�P E�� E�0 E� E�� FH� Fz  F�� F� Gp GCP E"� BM_�n+I���?6���H�?e�d����H �@BQ�,�� �JL˨?�>m:`Fn�f?m�<�=��F�~3BK��?�u�F�a ?�/=��G"�M?o���yGj�3BJ�]@C/GSN?��;��F_6�?l���JH
_�BJ�?�(�G�1�?˽�6�G�z ?k�A��h+H�ހBI>(����I��?KQ���H@?s�q�hH�� BL!H����J
�(?H�>�F��f?y5��z
GV/3BK������FW� ?	R�>K%�Gt�?w�Ͻ���G:�fBJ��?�W�GlVf?
�>/�G�� ?w佟[G��BJ�
@�G�π?	�=��pHe@?v�Q���G�m�BJ��@A�H�@?	��=�OH</�?w<����sH��BJ�b@Q�H� @?;�>l>uG6�3?v���3jF��fBJiD@L TI:j`?�7>6�OF��?uZ:����GV�BJ+Q@L �Ia��?�'>1W�F�A ?u����G% BJF�@�e�I_!@?*�>i>F�Gf?t�½�_=F��3BJ��@�y}Ig��?.>U��F�x�?q�׽��F��BJh�@}�I�Y�?��=ˑ F���?nsŽ�0�G�MBKbh�*�F[��?$�u>ٌF� ?p�1�5a�E� BJ?c@19I���?T�>��E��?bĽ}u\F1)BJ/ ?�EI�@�?1��>|��F(4 ?k+$����F_6�BJ�@�^J��?[S�>�:E�H?\-���%zFw��BJh@�.�I�� ?bn"=�t�E��R?X?սl|EFl� BJ	�@Z��Jo�?|<�<c��FQ� ?N�����G��BJ"�?�0Iɤ�?s�;8F_6�?L{���5HPv BJ0;?�N�I�#?�M9�RF_6�?K�}��:H�5�BI��?���I|5�?�m7��OF_6�?J���z�H�D�BCx�@(�G��?��5X�
F_6�?I����RH{ �BJ8R?�gI�J ?��5"l�F_6�?G���ˏH��`BJ*e?�[I�Ӑ?��4��0F_6�?D ����H�J�BJ��@!�.H�� ?�>7l�Gd�?u�׽���G 