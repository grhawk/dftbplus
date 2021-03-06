!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2017  DFTB+ developers group                                                      !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'

!!* Contains the HSD (Human readable Structured Data) parser.
!!* @desc
!!* <p>
!!*   The HSD format is a more or less user friendly input format, which can
!!*   be easily converted to a simplified XML format. The parser returns a
!!*   DOM-tree, which can be further processed. The returned tree contains
!!*   also information about the original name and position of the keywords
!!*   in the original HSD format, in order to enable user friendly error
!!*   messages, if inconsistent data is detected during the processing of the
!!*   DOM-tree.
!!* </p><p>
!!*   For the specification of the HSD format see the sample input
!!* </p>
module hsdparser
  use assert
  use message
  use fileid
  use charmanip
  use xmlutils
  use xmlf90
  implicit none
  private

  !!* Wrapper around the parsing function
  interface parseHSD
    module procedure parseHSD_stdin
    module procedure parseHSD_file
    module procedure parseHSD_opened
  end interface


  !!* Wrapper around the HSD dumping
  interface dumpHSD
    module procedure dumpHSD_file
    module procedure dumpHSD_opened
  end interface


  !! Main token separator characters
  integer, parameter :: nSeparator = 7
  character(len=*), parameter :: sIncludeXML = "<<!"
  character(len=*), parameter :: sIncludeParsed = "<<+"
  character(len=*), parameter :: sIncludeUnparsed = "<<<"
  character(len=*), parameter :: sSingleOpen = "=  "
  character(len=*), parameter :: sOpen = "{  "
  character(len=*), parameter :: sClose = "}  "
  character(len=*), parameter :: sSingleClose = ";  "
  character(len=*), parameter :: separators(nSeparator) = &
      &(/ sIncludeXML, sIncludeParsed, sIncludeUnparsed, &
      &sSingleOpen, sOpen, sClose, sSingleClose /)

  !! Other parsed characters
  character(len=*), parameter :: sModifierOpen = "["
  character(len=*), parameter :: sModifierClose = "]"
  character(len=*), parameter :: sComment = "#"

  !! Extension related stuff
  integer, parameter :: nExtension = 5
  character(len=*), parameter :: sExtendIfPresentOrDie = "+"
  character(len=*), parameter :: sExtendIfPresent = "?"
  character(len=*), parameter :: sExtendIfPresentOrCreate = "*"
  character(len=*), parameter :: sCreateIfNotPresent = "/"
  character(len=*), parameter :: sReplaceIfPresentOrCreate = "!"
  character(len=*), parameter :: extensions(nExtension) = &
      &(/ sExtendIfPresentOrDie, sExtendIfPresent, sExtendIfPresentOrCreate, &
      &sCreateIfNotPresent, sReplaceIfPresentOrCreate /)


  !! Name and file descriptor from standard input/output
  character(len=*), parameter :: stdin = "*"
  character(len=*), parameter :: stdout = "*"
  integer, parameter :: fdStdin = 0
  integer, parameter :: fdStdout = 0

  !! Forbidden (even if quoted) characters in the iput
  character(len=*), parameter :: forbiddenChars = "<>"

  !! Attribute names
  character(len=*), parameter :: attrStart = "start"
  character(len=*), parameter :: attrEnd = "end"
  character(len=*), parameter :: attrFile = "file"
  character(len=*), parameter :: attrName = "name"
  character(len=*), parameter :: attrModifier = "m"
  character(len=*), parameter :: attrList = "l"

  !! Length of a parsed line
  integer, parameter :: lc = 1024

  !! Maximal record lenght on output in characters (bytes).
  !! If text nodes bigger than that occur runtime error can be expected.
  integer, parameter :: MAXRECL = 1024 * 1024

  !! Name of the root tag
  character(len=lc) :: rootName

  !! Pointer to the top of the currently processed document
  !! (modified only in parseHSD and replaceTreeFromFile)
  type(fnode), pointer :: myDoc

  !! Format of the input line
  character(len=lc) :: lineFormat = ""


  public :: parseHSD, dumpHSD, dumpHSDAsXML, stdin, stdout, newline
  public :: getNodeHSDName, getHSDPath
  public :: attrStart, attrEnd, attrFile, attrName, attrModifier, attrList


contains

  !!* Parses HSD format from stdandard input
  !!* @param initRootName Name of the root tag of the resulting XML-tree
  !!* @param xmlDoc       DOM-tree of the parsed input on exit
  subroutine parseHSD_stdin(initRootName, xmlDoc)
    character(len=*), intent(in) :: initRootName
    type(fnode), pointer :: xmlDoc

    call parseHSD_opened(initRootName, fdStdin, stdin, xmlDoc)

  end subroutine parseHSD_stdin



  !!* Parser HSD format from a file
  !!* @param rootName Name of the root tag, which should contain the parsed tree
  !!* @param file     Name of the file (used in error messages)
  !!* @param xmlDoc       DOM-tree of the parsed input on exit
  subroutine parseHSD_file(initRootName, file, xmlDoc)
    character(len=*), intent(in) :: initRootName
    character(len=*), intent(in) :: file
    type(fnode), pointer :: xmlDoc

    integer, save :: fd = -1
    integer :: iostat

    if (fd == -1) then
      fd = getFileId()
    end if
    open(fd, file=file, iostat=iostat, status='old', action='read', recl=lc)
    if (iostat /= 0) then
      call parsingError("Error in opening file '" // trim(file) //"'.", &
          &file, -1)
    end if
    call parseHSD_opened(initRootName, fd, file, xmlDoc)
    close(fd, iostat=iostat)

  end subroutine parseHSD_file



  !!* Parses HSD format from an already opened file
  !!* @param rootName Name of the root tag, which should contain the parsed tree
  !!* @param fd       File descriptor of the open file containing the input
  !!* @param file     Name of the file (used in error messages)
  !!* @param xmlDoc       DOM-tree of the parsed input on exit
  subroutine parseHSD_opened(initRootName, fd, file, xmlDoc)
    character(len=*), intent(in) :: initRootName
    integer, intent(in) :: fd
    character(len=*), intent(in), optional :: file
    type(fnode), pointer :: xmlDoc

    type(fnode), pointer :: rootNode, dummy
    logical :: tFinished
    integer :: curLine
    character(len=lc) :: residual, curFile

    if (present(file)) then
      curFile = file
    elseif (fd == fdStdin) then
      curFile = stdin
    else
      curFile = "???"
    end if

    if (len_trim(lineFormat) == 0) then
      lineFormat = "(A" // i2c(lc) // ")"
    end if
    rootName = tolower(initRootName(:min(lc, len(initRootName))))
    myDoc => createDocumentNode()
    rootNode => createElement(trim(rootName))
    dummy => appendChild(myDoc, rootNode)
    curLine = 0
    residual = ""
    tFinished = parse_recursive(rootNode, 0, residual, .false., fd, curFile, &
        &0, curLine, &
        &(/ .true., .true., .true., .true., .true., .true., .true. /), .false.)
    xmlDoc => myDoc
    myDoc => null()

  end subroutine parseHSD_opened



  !!* Recursive parsing function for the HSD parser making the actual work
  !!* @param curNode     Node which should contain parsed input
  !!* @param depth       Number of open blocks/assignments.
  !!* @param residual    Unparsed text from the previous line
  !!* @param tRightValue Is next parsed token a right value of an assignment?
  !!* @param fd          File descriptor of the input
  !!* @param curFile     Name of the current input file
  !!* @param fileDepth   Number of open files
  !!* @param curLine     Number of current line in the current file
  !!* @param parsedTypes True for those separators, which should be parsed
  !!* @return            True, if parsing is done
  recursive function parse_recursive(curNode, depth, residual, tRightValue, &
      &fd, curFile, fileDepth, curLine, parsedTypes, tNew) result (tFinished)
    type(fnode), pointer             :: curNode
    integer, intent(in)              :: depth
    character(len=lc), intent(inout) :: residual
    logical, intent(in)              :: tRightValue
    integer, intent(in)              :: fd
    character(len=lc), intent(in)    :: curFile
    integer, intent(in)              :: fileDepth
    integer, intent(inout)           :: curLine
    logical, intent(in)              :: parsedTypes(nSeparator)
    logical, intent(in)              :: tNew
    logical                          :: tFinished

    character(len=lc) :: strLine, word

    type(fnode), pointer :: childNode, dummy
    type(string) :: buffer
    integer :: newFile
    integer :: iostat
    integer :: iType, sepPos
    integer :: newCurLine
    logical :: tTagClosed, tNewNodeCreated
    logical :: newParsedTypes(nSeparator)
    integer :: nTextLine
    integer :: iTmp
    integer :: nodetype

    tTagClosed = .false.
    tFinished = .false.
    tNewNodeCreated = .false.
    newFile = -1
    nTextLine = 0
    nodetype = 0

    lpMain: do while ((.not. tTagClosed) .and. (.not. tFinished))

      !! Read in next line or process residual from last line.
      if (len_trim(residual) /= 0) then
        strLine = adjustl(residual)
        residual = ""
      else
        if (fd == fdStdin) then
          read (*, trim(lineFormat), iostat=iostat) strLine
        else
          read (fd, trim(lineFormat), iostat=iostat) strLine
        end if
        curLine = curLine + 1
        call convertWhitespaces(strLine)
        strLine = adjustl(strLine)
        !! If reading error (e.g. EOF) -> close current scope
        if (iostat /= 0) then
          tTagClosed = .true.
          if (depth /= 0) then
            call getAttribute(curNode, attrStart, buffer)
            call parsingError("Unexpected end of input (probably open node at &
                &line " // char(buffer) // " or after).", curFile, curLine)
          end if
          !! If outermost file, we are ready
          if (fileDepth == 0) then
            tFinished = .true.
          end if
          exit
        end if
      end if

      !! Remove comments
      iTmp = unquotedIndex(strLine, sComment)
      if (iTmp /= 0) then
        strLine = strLine(:iTmp-1)
      end if

      !! Get first occurance of any separator
      call getFirstOccurance(strLine, separators, parsedTypes, iType, sepPos)

      !! Handle various closing operators.
      select case (iType)
      case (6)
        !! Block closing char on level zero is invalid
        if (depth == 0) then
          call parsingError("Invalid block closing sign.", curFile, curLine)
        end if
        !! If block closing char is not first char of the line, text before it
        !! will be appended as text, and residual line reparsed in next cycle
        if (sepPos == 1) then
          tTagClosed = .true.
          residual = strLine(sepPos+1:)
        else
          iType = 0
          residual = strLine(sepPos:)
          strLine = strLine(:sepPos-1)
        end if
      case (7)
        if (tRightValue) then
          !! Single assignment is terminated by a semi-colon. Text after
          !! semicolon must be reparsed in next cycle.
          iType = 0
          residual = strLine(sepPos+1:)
          strLine = strLine(:sepPos-1)
          tTagClosed = .true.
        else
          call parsingError("Invalid assignment separator", curFile, curLine)
        end if
      end select

      !! Ignore empty lines
      if (len_trim(strLine) == 0) then
        cycle lpMain
      end if

      !! Check for forbidden characters in the current line
      if (.not. (iType == 1 .or. iType == 2 .or. iType == 3)) then
        call checkForbiddenChars(strLine, curFile, curLine)
      end if

      !! Process non-closing separators
      select case (iType)
      case(0)
        if (nodetype > 0) then
          call parsingError("Node already contains subnodes, no text content&
              & allowed any more", curFile, curLine)
        end if
        !! No separator found -> Add entire line as text
        !! If current node already contains text, prepend newline before
        !! appending. (Otherwise newlines in the text would get lost.)
        if (associated(curNode)) then
          nTextLine = nTextLine + 1
          if (nTextLine > 1) then
            childNode => createTextNode(newline // trim(strLine))
          else
            childNode => createTextNode(trim(strLine))
          end if
          dummy => appendChild(curNode, childNode)
        end if
        nodetype = -1

      case(1)
        !! XML inclusion
        if (associated(curNode)) then
          if (sepPos /= 1) then
            call parsingError("Invalid character before file inclusion &
                &operator", curFile, curLine)
          end if
          strLine = adjustl(unquote(strLine))
          word = adjustl(strLine(len(sIncludeXML)+1:len_trim(strLine)))
          if (len_trim(word) /= 0) then
            if (depth == 0 .and. fileDepth == 0 &
                &.and. .not. associated(getFirstChild(curNode))) then
              call replaceTreeFromFile(curNode, trim(word))
            else
              call parsingError("XML inclusion must be the first statement &
                  &in the first input file.", curFile, curLine)
            end if
          else
            call parsingError("No file name specified after the inclusion &
                &operator.", curFile, curLine)
          end if
        end if


      case(2, 3)
        !! File inclusion operator -> append content of new file to current node
        if (associated(curNode)) then
          if (sepPos /= 1) then
            call parsingError("Invalid character before file inclusion &
                &operator", curFile, curLine)
          end if
          strLine = adjustl(unquote(strLine))
          if (iType == 2) then
            word = adjustl(strLine(len(sIncludeParsed)+1:len_trim(strLine)))
          else
            word = adjustl(strLine(len(sIncludeUnparsed)+1:len_trim(strLine)))
          end if

          if (len_trim(word) == 0) then
            call parsingError("No file name specified after the inclusion &
                &operator.", curFile, curLine)
          end if

          if (newFile == -1) then
            newFile = getFileID()
          end if
          open(newFile, file=trim(word), status='old', action='read', &
              &iostat=iostat, recl=lc)
          if (iostat /= 0) then
            call parsingError("Error in opening file '" // trim(word) // &
                &"'.", curFile, curLine)
          end if
          strLine = ""
          newCurLine = 0
          if (iType == 2) then
            !! Everything is parsed
            newParsedTypes = (/ .true., .true., .true., .true., .true., &
                &.true., .true. /)
          else
            !! Nothing is parsed
            newParsedTypes = (/ .false., .false., .false., .false., .false., &
                &.false., .false. /)
          end if
          tFinished = parse_recursive(curNode, 0, strLine, .false., newFile, &
              &word, fileDepth + 1, newCurLine, newParsedTypes, .false.)
          close(newFile, iostat=iostat)
        end if

      case(4)
        !! Assignment
        if (nodetype < 0) then
          call parsingError("Node already contains free text, no child nodes&
              & allowed any more", curFile, curLine)
        end if
        word = adjustl(strLine(:sepPos-1))
        strLine = adjustl(strLine(sepPos+1:))
        if (len_trim(word) == 0) then
          call parsingError("Missing field name on the left of the &
              &assignment.", curFile, curLine)
        elseif (len_trim(strLine) == 0) then
          call parsingError("Missing value on the right of the assignment.", &
              &curFile, curLine)
        elseif (nTextLine > 0) then
          call parsingError("Unparsed text before current node", curFile, &
              &curLine)
        end if
        if (associated(curNode)) then
          childNode => createChildNode(curNode, word, curLine, curFile)
          tNewNodeCreated = .true.
        else
          childNode => null()
        end if
        !! Only block opening/closing sign and single child separator are parsed
        newParsedTypes = (/ .false., .false., .false., .false., .true., &
            &.true., .true. /)
        tFinished = parse_recursive(childNode, depth+1, strLine, .true., fd, &
            &curFile, fileDepth, curLine, newParsedTypes, tNewNodeCreated)
        residual = strLine
        nodetype = 1


      case(5)
        if (nodetype < 0) then
          call parsingError("Node already contains free text, no child nodes&
              & allowed any more", curFile, curLine)
        end if
        !! Block opening sign
        word = adjustl(strLine(:sepPos-1))
        strLine = adjustl(strLine(sepPos+1:))
        if (nTextLine > 0) then
          call parsingError("Unparsed text before current node", curFile, &
              &curLine)
        end if
        if (associated(curNode)) then
          ! Currently node without name is allowed to support "= {" construct
          ! Should be turned to parsing error to deprecate that construct.
          if (len_trim(word) == 0) then
            childNode => curNode
            call setAttribute(curNode, attrList, "")
            !call parsingError("Node without name not allowed.", curFile,&
            !    & curLine)
          else
            childNode => createChildNode(curNode, word, curLine, curFile)
            tNewNodeCreated = .true.
          end if
        else
          childNode => null()
        end if
        newParsedTypes = (/ .true., .true., .true., .true., .true., .true., &
            &.true. /)
        tFinished = parse_recursive(childNode, depth+1, strLine, .false., &
            &fd, curFile, fileDepth, curLine, newParsedTypes, tNewNodeCreated)
        residual = strLine
        nodetype = 1

      end select

      tTagClosed = tTagClosed .or. tRightValue

    end do lpMain

    !! Set end attribute on tag end and normalise text nodes
    if (tTagClosed .and. associated(curNode)) then
      if (tNew) then
        call setAttribute(curNode, attrEnd, i2c(curLine))
      end if
      if (nTextLine > 1) then
        call normalize(curNode)
      end if
    end if

  end function parse_recursive



  !!* Creates a child node with attributes related to the HSD input.
  !!* @param parentNode Parent node containing of the child to be created
  !!* @param childName  Name of the new child
  !!* @param curLine    Number of the current line
  !!* @param file       Name of the current file
  !!* @return Pointer to the new (appended) child node
  function createChildNode(parentNode, childName, curLine, file) &
      &result(newChild)
    type(fnode), pointer          :: parentNode
    character(len=lc), intent(in) :: childName
    integer, intent(in)           :: curLine
    character(len=lc), intent(in) :: file
    type(fnode), pointer          :: newChild

    type(fnode), pointer :: dummy, sameChild
    character(len=lc) :: lowerName, truncName, modifier
    logical :: tModifier, tCreate
    integer :: pos1, pos2, iType
    integer :: ii

    truncName = childName
    lowerName = tolower(childName)

    !! Look for any extension operator
    iType = 0
    pos1 = 0
    do ii = 1, nExtension
      pos1 = len_trim(extensions(ii))
      if (lowerName(:pos1) == trim(extensions(ii))) then
        iType = ii
        exit
      end if
    end do

    !! Cut extension operator from the field name
    if (iType /= 0) then
      lowerName = lowerName(pos1+1:)
      truncName = truncName(pos1+1:)
    end if

    !! Look for modifier after field name
    tModifier = .false.
    pos1 = index(lowerName, sModifierOpen)
    pos2 = index(lowerName, sModifierClose)
    if (pos1 == 0) then
      if (pos2 /= 0) then
        call parsingError("Unbalanced modifier opening sign.", file, curLine)
      end if
    else
      if (pos2 /= len_trim(lowerName)) then
        call parsingError("Invalid character(s) after modifier closing sign.",&
            &file, curLine)
      end if
      !! Remove modifier from field name
      modifier = adjustl(truncName(pos1+1:pos2-1))
      lowerName = adjustl(lowerName(:pos1-1))
      truncName = adjustl(truncName(:pos1-1))
      tModifier = .true.
    end if

    !! Check if field name is nonempty
    if (len_trim(lowerName) == 0) then
      call parsingError("Missing field name", file, curLine)
    end if

    !! Create child according extension operator
    tCreate = .false.
    if (iType == 0) then
      !! Create and append new node
      tCreate = .true.
    else
      !! Look for already present node with the same name
      sameChild => getFirstChildByName(parentNode, trim(lowerName))
      if (associated(sameChild)) then
        !! We have found a block with the same name
        if (iType == 4) then
          newChild => null()
          return
        elseif (iType == 5) then
          dummy => removeChild(parentNode, sameChild)
          call destroyNode(sameChild)
          tCreate = .true.
        else
          newChild => sameChild
        end if
      else
        !! We did not found a child with the same name
        select case (iType)
        case(1)
          call parsingError("Containing block does not contain a(n) '" &
              &// trim(truncName) // "' block yet.", file, curLine)
        case(2)
          newChild => null()
          return
        case(3,4,5)
          tCreate = .true.
        end select
      end if
    end if

    !! Create and append the node
    if (tCreate) then
      newChild => createElement(trim(lowerName))
      dummy => appendChild(parentNode, newChild)
    end if

    !! Set useful attributes
    call setAttribute(newChild, attrStart, i2c(curLine))
    call setAttribute(newChild, attrName, trim(truncName))
    call setAttribute(newChild, attrFile, trim(file))
    if (tModifier) then
      call setAttribute(newChild, attrModifier, trim(modifier))
    end if

  end function createChildNode



  !!* Checks for forbidden characters and issue error message, if any found.
  !!* @param str     String to investigate
  !!* @param curFile Name of the current file
  !!* @param curLine Number of the current line in the current file
  subroutine checkForbiddenChars(str, curFile, curLine)
    character(len=*), intent(in) :: str
    character(len=*), intent(in) :: curFile
    integer, intent(in) :: curLine

    if (scan(str, forbiddenChars) /= 0) then
      call parsingError("Invalid character(s).", curFile, curLine)
    end if

  end subroutine checkForbiddenChars



  !!* Issues a parsing error message containing file name and line number.
  !!* @param message Parsing error message
  !!* @param file    Name of the current file
  !!* @param line    Number of current line
  subroutine parsingError(message, file, line)
    character(len=*), intent(in) :: message
    character(len=lc), intent(in) :: file
    integer, intent(in) :: line

    character(len=lc) :: msgArray(2)

    if (trim(file) == stdin) then
      write (msgArray(1), 9990) line
9990  format("HSD parser error: Standard input, Line ",I5,".")
    else
      !! Watch out to trunk away enough from the file name to prevent overflow
      write (msgArray(1), 9991) trim(file(1:lc-40)), line
9991  format("HSD parser error: File '",A,"', Line",I5,".")
    end if
    write (msgArray(2), "(A)") trim(message(:min(lc, len(message))))
    call error(msgArray)

  end subroutine parsingError



  !!* Dumps the DOM-tree of a HSD document to a file.
  !!* @param myDoc    DOM-tree of a HSD document
  !!* @param fileName File for the XML-dump.
  !!* @descr This routine pretty prints the XML-tree in the specified file.
  !!*   Attributes related to the HSD document (e.g. line number, file etc.)
  !!*   are not printed.
  subroutine dumpHSDAsXML(myDoc, fileName)
    type(fnode), pointer :: myDoc
    character(len=*), intent(in) :: fileName

    type(xmlf_t) :: xf
    type(fnode), pointer :: fp

    call xml_OpenFile(fileName, xf, indent=.true.)
    call xml_AddXMLDeclaration(xf)
    fp => getFirstChild(myDoc)
    if (associated(fp)) then
      call dumpHSDAsXML_recursive(xf, fp)
    end if
    call xml_Close(xf)

  end subroutine dumpHSDAsXML



  !!* Recursive working horse for dumpHSD
  !!* @param xf   XML pretty printer data
  !!* @param node Node to prety print
  recursive subroutine dumpHSDAsXML_recursive(xf, node)
    type(xmlf_t), intent(inout) :: xf
    type(fnode), pointer :: node

    type(string) :: txt, name, value
    type(fnode), pointer :: fp
    type(fnamedNodeMap), pointer :: attribs
    integer :: ii

    call getNodeName(node, txt)
    if (getNodeType(node) == TEXT_NODE) then
      call getNodeValue(node, value)
      call xml_AddPCData(xf, trim2(char(value)))
    else
      call xml_NewElement(xf, char(txt))
      attribs => getAttributes(node)
      do ii = 0, getLength(attribs) - 1
        fp => item(attribs, ii)
        call getNodeName(fp, name)
        !! Currently only the modifier and single attributes are dumped
        if (name == attrModifier .or. name == attrList) then
          call getNodeValue(fp, value)
          call xml_AddAttribute(xf, char(name), char(value))
        end if
      end do
      fp => getFirstChild(node)
      do while (associated(fp))
        call dumpHSDAsXML_recursive(xf, fp)
        fp => getNextSibling(fp)
      end do
      call xml_EndElement(xf, char(txt))
    end if

  end subroutine dumpHSDAsXML_recursive



  !!* Replaces the tree
  !!* @param curNode Node, which should contained the parsed children from file
  !!* @param file    File to parse
  !!* @note The solution with access to a global module variable is not very
  !!*   elegant, but it saves the deep cloning of the parsed document.
  subroutine replaceTreeFromFile(curNode, file)
    type(fnode), pointer :: curNode
    character(len=*), intent(in) :: file

    type(fnode), pointer :: newDoc, rootNode

    newDoc => parsefile(file)
    call removeSpace(newDoc)
    call normalize(newDoc)
    rootNode => getLastChildByName(newDoc, trim(rootName))
    if (.not. associated(rootNode)) then
      call parsingError("File '" // file // "' does not contain '" // &
          & trim(rootName) // "' node.", file, -1)
    else
      call destroyNode(myDoc)
      myDoc => newDoc
      curNode => rootNode
    end if

  end subroutine replaceTreeFromFile



  !!* Dumps a HSD tree in a file.
  !!* @param myDoc The DOM tree
  !!* @param file  Name of the file
  subroutine dumpHSD_file(myDoc, file)
    type(fnode), pointer :: myDoc
    character(len=*), intent(in) :: file

    integer, save :: fd = -1
    integer :: iostat
    character(len=lc) :: fileName

    if (fd == -1) then
      fd = getFileId()
    end if
    open(fd, file=file, iostat=iostat, status='replace', action='write',&
        &recl=MAXRECL)
    if (iostat /= 0) then
      fileName = file
      call parsingError("Error in opening file for the HSD output.", fileName, &
          &-1)
    end if
    call dumpHSD_opened(myDoc, fd)
    close(fd)

  end subroutine dumpHSD_file



  !!* Dumps a DOM-tree representing a HSD input in HSD format to an opened file.
  !!* @param myDoc The DOM tree
  !!* @param fd    File descriptor for an open file where output should go.
  subroutine dumpHSD_opened(myDoc, fd)
    type(fnode), pointer :: myDoc
    integer, intent(in) :: fd

    type(fnode), pointer :: rootNode
    type(fnode), pointer :: child
    type(string) :: buffer

    rootNode => getFirstChild(myDoc)
    if (.not. associated(rootNode)) then
      return
    end if
    child => getFirstChild(rootNode)
    do while (associated(child))
      call dumpHSD_recursive(child, 0, fd, .false., buffer)
      child => getNextSibling(child)
    end do

  end subroutine dumpHSD_opened



  !!* Recursive working horse for the dumpHSD routine.
  !!* @param node        Node to dump
  !!* @param indent      Current indentation level
  !!* @param fd          File descriptor for an open file where output should go
  !!* @param tRightValue Is current node the right hand side of an assignment?
  !!* @param buffer      Buffer for storing temporary strings
  recursive subroutine dumpHSD_recursive(node, indent, fd, tRightValue, buffer)
    type(fnode),      pointer       :: node
    integer,          intent(in)    :: indent
    integer,          intent(in)    :: fd
    logical,          intent(in)    :: tRightValue
    type(string),     intent(inout) :: buffer

    type(fnode), pointer :: child, attr
    logical :: tOpenBlock

    !! Text nodes are printed including trailing newline. No further processing
    if (getNodeType(node) == TEXT_NODE) then
      call getNodeValue(node, buffer)
      write (fd, "(A)", advance="yes") trim2(char(buffer))
      return
    end if

    !! If left value -> indent
    if (.not. tRightValue) then
      write (fd, "(A)", advance="no") repeat(" ", indent)
    end if

    !! Get and print tag name
    attr => getAttributeNode(node, attrName)
    if (associated(attr)) then
      call getNodeValue(attr, buffer)
    else
      call getNodeName(node, buffer)
    end if
    write (fd, "(A)", advance="no") char(buffer)

    !! Get and print modifier
    attr => getAttributeNode(node, attrModifier)
    if (associated(attr)) then
      call getNodeValue(attr, buffer)
      if (len(buffer) > 0) then
        write (fd, "(' ',A,A,A)", advance="no") trim(sModifierOpen), &
            &char(buffer), trim(sModifierClose)
      end if
    end if

    child => getFirstChild(node)
    if (tRightValue) then
      if (associated(child)) then
        !! Children exist -> open block
        write (fd, "(' ',A)", advance="yes") trim(sOpen)
        tOpenBlock = .true.
      else
        !! No children -> Open and close block immediately and return
        write (fd, "(' ',A,A)", advance="yes") trim(sOpen), trim(sClose)
        return
      end if
    else
      !! We are the left hand side of an assignment -> write assignment sign
      write (fd, "(' ',A,' ')", advance="no") trim(sSingleOpen)
      attr => getAttributeNode(node, attrList)
      if (associated(child)) then
        tOpenBlock = .false.
        if (associated(attr) .or. associated(getNextSibling(child))) then
          !! RHS has many children or is signed as block -> open a block
          tOpenBlock = .true.
        elseif (getNodeType(child) == TEXT_NODE) then
          !! RHS's child is text -> open block if it contains several tokens
          call getNodeValue(child, buffer)
          tOpenBlock = (unquotedScan(trim2(char(buffer)), whiteSpaces) /= 0)
        end if
        if (tOpenBlock) then
          write (fd, "(A)", advance="yes") trim(sOpen)
        end if
      else
        !! No child, write an empty block and return
        write (fd, "(A,A)", advance="yes") trim(sOpen), trim(sClose)
        return
      end if
    end if

    !! Process children
    do while (associated(child))
      if (tOpenBlock) then
        call dumpHSD_recursive(child, indent+2, fd, .false., buffer)
      else
        call dumpHSD_recursive(child, indent, fd, .true., buffer)
      end if
      child => getNextSibling(child)
    end do

    !! Dump closing sign
    if (tOpenBlock) then
      write (fd, "(A)", advance="no") repeat(" ", indent)
      write (fd, "(A)", advance="yes") trim(sClose)
    end if

  end subroutine dumpHSD_recursive



  !!* Returns the name of a node, if present in pretty printing format.
  !!* @param node Node to investigate.
  subroutine getNodeHSDName(node, name)
    type(fnode), pointer :: node
    type(string), intent(inout) :: name

    call getAttribute(node, attrName, name)
    if (len(name) == 0) then
      call getNodeName(node, name)
    end if

  end subroutine getNodeHSDName



  !!* Returns the path of a node, if possible in pretty printing format.
  !!* @param node        Node to investigate
  !!* @param path        String containing the path on return
  !!* @param excludeRoot If root node should be excluded
  subroutine getHSDPath(node, path, excludeRoot)
    type(fnode), pointer :: node
    type(string), intent(inout) :: path
    logical, intent(in), optional :: excludeRoot

    type(fnode), pointer :: parent
    type(string) :: buffer
    logical :: inclRoot

    if (present(excludeRoot)) then
      inclRoot = .not. excludeRoot
    else
      inclRoot = .true.
    end if

    call getNodeHSDName(node, path)
    parent => getParentNode(node)
    if (associated(parent)) then
      call getHSDPath_recursive(parent, path, inclRoot, buffer)
    elseif (.not. inclRoot) then
      path = ""
    end if

  end subroutine getHSDPath



  !!* Working horse for the getHSDPath routine
  !!* @param node     Node to look for
  !!* @param path     String containing the path until now
  !!* @param inclRoot If document root should be included
  !!* @param buffer   Buffer for strings (to avoid destruction at every call)
  recursive subroutine getHSDPath_recursive(node, path, inclRoot, buffer)
    type(fnode), pointer :: node
    type(string), intent(inout) :: path
    logical, intent(in) :: inclRoot
    type(string), intent(inout) :: buffer

    type(fnode), pointer :: parent

    parent => getParentNode(node)
    if (associated(parent) .or. inclRoot) then
      call prepend_to_string(path, "/")
      call getNodeHSDName(node, buffer)
      call prepend_to_string(path, buffer)
      if (associated(parent)) then
        call getHSDPath_recursive(parent, path, inclRoot, buffer)
      end if
    end if

  end subroutine getHSDPath_recursive



end module hsdparser
