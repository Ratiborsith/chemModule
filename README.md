# Программный модуль анализа химических формул в патентных документах

## Цель
Цель настоящей работы заключается в разработке программного модуля, способного распознавать и классифицировать информацию, а затем преобразовывать полученные данные в единый формат Fingerprint для последующего сравнения. Настоящая работа направлена на создание программного продукта, способного автоматизировать процессы анализа и сравнения формул из патентных документов. 


## Задачи
Для достижения данной цели необходимо решить следующие задачи:
1. Провести анализ предметной области, рассмотреть системы-аналоги, методы решения задачи распознавания графических и текстовых представлений химических формул, преобразования различных типов формул (IUPAC, SMILES, Inchi, MOL, структурная формула).
2. Разработать алгоритм распознавания графических представлений химических формул.
3. Разработать алгоритм распознавания текстовых представлений химических формул.
4. Разработать алгоритм преобразований различных представлений химических формул.
5. Спроектировать программный модуль.
6. Программно реализовать модуль и проверить его эффективность.

## Установка
Чтобы запустить этот программный модуль анализа химических формул в патентных документах, убедитесь, что на вашей машине(MacOs Sonoma или Windows 10+) установлен Python 3.10 или новее. Клонируйте этот репозиторий, а затем установите необходимые зависимости:
```bash
git clone https://github.com/Ratiborsith/chemModule
```

## Использование
Для того, чтобы запустить графический интерфейс программного модуля введите следующую команду в терминале:
```bash
python app.py
```
После запуска вас встретит простой и интуитивно понятный интерфейс. Вы сможете осуществлять следующие функции:
1.	 Сравнение по коэффициенту Tanimoto химических формул разных представлений друг с другом. Доступны для сравнения следующие химические представления: InChi, MOLFILE, SMILES, структурная формула, IUPAC
2.	Сравнение патентных документов. Сравниваются по парам все химические формулы, которые присутствуют в патенте.
3.	Парсинг структурных и IUPAC представлений из Google Patents
4.	Визуализация представленной информации в удобном для пользователя формате. 


## Технологии
#### • 	СУБД: SQLite3 для хранения данных.
#### •	Библиотеки для химического анализа: RDKit, Cirpy.
#### •	Веб-сервер: Flask.
#### •	Средства аутентификации и авторизации: Flask-Login.
#### •	Скрипты интерфейса: JavaScript, ajax
#### •	Внешний вид веб-интерфейса: html, css, scss. Для работы с формами: flask wtforms.
#### •	Библиотеки для подмодуля парсинга и краулинга: Selenium, BeautifulSoup.
#### •	Используемая модель для конвертации структурных формул: OSRA.

## Архитектура
<img width="397" alt="image" src="https://github.com/Ratiborsith/chemModule/assets/85187788/2aceb93e-ce1f-4e14-b350-1a5c366853e4"> <br>
1. – обращение к базе Google Patents по ссылке патента;
2. – получение из Google Patents кода страницы патента;
3. – вызов блока парсинга формул;
4. – передача в БД спарсенных химических представлений IUPAC и структурной формулы;
5. – получение из БД списка патентов;
6. – получение из БД списка патентов и содержащихся в каждом выбранном патенте химических соединений;
7. – вызов блока создания молекул RDKit;
8. – вызов блока получения fingerprint;
9. – Обращение к подмодулю сравнения. Вызов блока сравнения формул по коэффициенту Tanimoto для двух fingerprint;
10. – вызов блока сравнения патентов;
11. – обращение к подмодулю визуализации. вызов блока результатов рассчета коэффициента Tanimoto;
12. – обращение к БД для проверки авторизации;
13. – получение из БД пароля для указанного по логину пользователя. Сравнение через Flask-Login с указанными;
14. – получение из БД списка патентов.




