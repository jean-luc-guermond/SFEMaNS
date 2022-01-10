// ---------------------------------------------------------------------
// $Id$
//
//    This file is part of SFEMaNS.
//
//    SFEMaNS is free software: you can redistribute it and/or modify
//    it under the terms of the GNU General Public License as published by
//    the Free Software Foundation, either version 3 of the License, or
//    (at your option) any later version.
//
//    SFEMaNS is distributed in the hope that it will be useful,
//    but WITHOUT ANY WARRANTY; without even the implied warranty of
//    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//    GNU General Public License for more details.
//
//    You should have received a copy of the GNU General Public License
//    along with SFEMaNS.  If not, see <http://www.gnu.org/licenses/>.
// ---------------------------------------------------------------------



/**
 * @page doc_SFEMaNS_read_user_data Fortran file read_user_data.f90

The file <tt>read_user_data.f90</tt> allows to add new datas in the code SFEMaNS. These datas are stocked in the derived data type called user. As the datas defined in the code SFEMaNS, these user datas are read in the <tt>data</tt> file described in this <a href='doc_SFEMaNS_data.html'><code>section</code></a>. A template is provided in the following directory: <tt>($SFEMaNS_DIR)/TEMPLATE</tt>.


To add new data, that will be read by the code SFEMaNS in the <tt>data</tt> file, you need to edit the file <tt>read_user_data.f90</tt> as follows.
<ol>
<li>The new data has to be declared in the module user_data_module. For example to define a real number called my_real, the following line is added.
\code
     REAL(KIND=8)           :: my_real
\endcode
We note that the user should only modify the module user_data_module between the following code lines:
\code
     !===I declare my own data here==================================================
     LOGICAL                                 :: if_my_stuff
     !.......Continue here ................................
     !===End I declare my own data here==============================================
\endcode
<li>The new data needs to be read by the code SFEMaNS. This is done in the subroutine <tt>read_user_data</tt> of the module user_data. To do so, we use the function <tt>find_string</tt> as follows:
\code
    CALL find_string(unit_file, '===Value of my_real', test)
    IF (test) THEN
       READ (unit_file, *) user%my_real
    ELSE
       user%my_real=0.d0
    END IF
\endcode
The inputs are unit_file (that represents the <tt>data</tt> file) and the line that is read by the code SFEMaNS in the <tt>data</tt> file (here it is '===Value of my_real'). The ouput is a logical denoted test. It is true if the line '===Value of my_real' is present in the <tt>data</tt> file and false otherwise.

We note that the user should only modify the subroutine <tt>read_user_data</tt> between the following code lines:
\code
    !===I add lines that the code SFEMaNS reads in the data file=========================
    !.......Continue here ................................
    !===End I add lines that the code SFEMaNS reads in the data file=====================
\endcode
</ol>

Remarks:
<ol>
<li>The real number my_real is stocked in user\%my_real.
<li>To have acces to the derived data type user in a function or a subroutine, you need to add the following module:
\code
   USE user_data
\endcode
</ol>






 */
